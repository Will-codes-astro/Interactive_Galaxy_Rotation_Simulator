"""The site's wording, outside the code, so Will can write it himself.

The text lives in one place only: HTML_TEMPLATE in generate_interactive_html.py.
This file reads every piece of visible text out of it (page HTML and the
sentences in the page's scripts, such as the guided tour cards) and can write
edited text back in. The generator calls write_snapshot() on every run, so
content/site_text.json is always the text of the site as generated. That
snapshot feeds the "Site text editor" artifact, where Will edits; his edits are
pulled back with `apply`. See .claude/skills/site-text/SKILL.md for the routine.

    python site_text.py apply EDITS.json   write edits into the generator
    python site_text.py check              Node syntax check of index.html's scripts
    python site_text.py status             is the editor artifact up to date?
    python site_text.py mark-published     record that the editor now shows this snapshot
    python site_text.py hook               Claude Code hook: after a commit, remind about the editor
"""

from __future__ import annotations

import hashlib
import html
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
GENERATOR_PATH = BASE_DIR / "generate_interactive_html.py"
INDEX_PATH = BASE_DIR / "index.html"
SNAPSHOT_PATH = BASE_DIR / "content" / "site_text.json"
PUBLISHED_PATH = BASE_DIR / "content" / "site_text_published.txt"
TEMPLATE_OPENING = 'HTML_TEMPLATE = r"""'

# Tags that sit inside a line of text; any other tag ends a piece of text.
INLINE_TAGS = {"a", "abbr", "b", "bdi", "br", "cite", "code", "em", "i", "kbd", "mark", "q",
               "s", "small", "span", "strong", "sub", "sup", "u", "var", "wbr", "time"}
VOID_TAGS = {"br", "wbr"}
HEADING_TAGS = {"h1", "h2", "h3", "h4", "h5", "h6"}
MIN_LETTERS_HTML = 2  # a piece of page text needs at least this many letters
MIN_WORDS_JS = 3  # a script string needs this many words to count as a sentence
# Script strings that are code, not sentences (CSS selectors, colours, units).
JS_CODE_HINTS = re.compile(r"^[.#\[]|:not\(|=>|;\S|;\s*$|[{}]|^rgba?\(|^\w+\(|\d+px")
# Script strings like 'curvePlot:Watch the cyan line' start with a code name the page reads.
CODE_PREFIX_RE = re.compile(r"^[a-z]\w*:(?! )")
MATH_TAGS = {"mi", "mo", "mn"}  # single maths symbols, not wording
ALLOWED_JS_ESCAPES = {"\\'", '\\"', "\\`"}

TAG_RE = re.compile(r"<!--.*?-->|<(/?)([a-zA-Z][\w-]*)\b[^>]*?(/?)>", re.S)
PAGE_RE = re.compile(r'data-page="([\w-]+)"')
JS_SECTION_RE = re.compile(r"\{id:'([\w-]+)'|function\s+(\w+)|const\s+([A-Z_]{3,})\s*=")
WORD_RE = re.compile(r"[A-Za-z]{2,}")
PLACEHOLDER_RE = re.compile(r"__[A-Z0-9_]+__")
EXPRESSION_RE = re.compile(r"\$\{[^}]*\}")


def read_template() -> tuple[str, int]:
    """The generator's source and the offset where HTML_TEMPLATE's text starts."""
    source = GENERATOR_PATH.read_bytes().decode("utf-8")
    start = source.index(TEMPLATE_OPENING) + len(TEMPLATE_OPENING)
    return source, start


def template_text(source: str, start: int) -> str:
    return source[start:source.index('"""', start)]


def plain(fragment: str) -> str:
    return html.unescape(re.sub(r"<[^>]+>", "", fragment))


# --- reading the page HTML ---------------------------------------------------

def _html_pieces(template: str):
    """Yield (start, end, tag, page) for each run of inline text between block tags,
    plus the (start, end) of every inline script, which _js_pieces reads."""
    scripts: list[tuple[int, int]] = []
    pieces = []
    run_start = None
    block_tag = "body"
    div_depth = 0
    page = "top bar"
    page_depth = None
    pos = 0
    for match in TAG_RE.finditer(template):
        if match.start() < pos:
            continue
        closing, name, _ = match.group(1), (match.group(2) or "").lower(), match.group(3)
        if match.group(0).startswith("<!--"):
            continue
        if name in INLINE_TAGS:
            if run_start is None:
                run_start = match.start()
            continue
        if run_start is not None:
            pieces.append((run_start, match.start(), block_tag, page))
            run_start = None
        if name in ("script", "style") and not closing:
            end = template.index(f"</{name}>", match.end())
            if name == "script" and "src=" not in match.group(0):
                scripts.append((match.end(), end))
            pos = end
            continue
        if name == "div":
            if closing:
                div_depth -= 1
                if page_depth is not None and div_depth < page_depth:
                    page, page_depth = "shared (menus, tour card, pop-ups)", None
            else:
                div_depth += 1
                found = PAGE_RE.search(match.group(0))
                if found and 'class="page' in match.group(0):
                    page, page_depth = found.group(1), div_depth
        if not closing:
            block_tag = name
        # Text right after this tag starts a new run.
        run_start = match.end()
    return pieces, scripts


def _trim(template: str, start: int, end: int) -> tuple[int, int] | None:
    """Shrink a run to whole, balanced text; None if it holds no real text."""
    fragment = template[start:end]
    left = len(fragment) - len(fragment.lstrip())
    right = len(fragment.rstrip())
    start, end = start + left, start + right
    # Drop closing tags at the front and opening tags at the end: they belong to outer text.
    changed = True
    while changed and start < end:
        changed = False
        head = re.match(r"\s*</[a-z]+>\s*", template[start:end])
        if head:
            start += head.end()
            changed = True
        tail = re.search(r"\s*<[a-z][^>]*>\s*$", template[start:end])
        if tail and not tail.group(0).strip().startswith("<br"):
            end = start + tail.start()
            changed = True
    fragment = template[start:end]
    if len(re.findall(r"[A-Za-z]", plain(fragment))) < MIN_LETTERS_HTML:
        return None
    opened = [t.lower() for t in re.findall(r"<([a-zA-Z]+)", fragment) if t.lower() not in VOID_TAGS]
    shut = [t.lower() for t in re.findall(r"</([a-zA-Z]+)", fragment)]
    if sorted(opened) != sorted(shut):
        return None
    return start, end


# --- reading the page scripts -------------------------------------------------

def _js_strings(script: str, offset: int):
    """Yield (start, end, quote) of every string literal in a script, offsets in the template."""
    i, n = 0, len(script)

    while i < n:
        c = script[i]
        if c == "/" and script.startswith("//", i):
            i = script.find("\n", i)
            i = n if i < 0 else i
        elif c == "/" and script.startswith("/*", i):
            i = script.index("*/", i) + 2
        elif c in "'\"`":
            j = i + 1
            while j < n and script[j] != c:
                if script[j] == "\\":
                    j += 1
                elif c == "`" and script.startswith("${", j):
                    depth, j = 1, j + 2
                    while depth:  # skip the expression, including strings inside it
                        if script[j] in "'\"`":
                            q, j = script[j], j + 1
                            while script[j] != q:
                                j += 2 if script[j] == "\\" else 1
                        depth += {"{": 1, "}": -1}.get(script[j], 0)
                        j += 1
                    continue
                elif c != "`" and script[j] == "\n":
                    break
                j += 1
            yield offset + i, offset + j + 1, c
            i = j + 1
        elif c == "/" and re.search(r"[(,=:\[!&|?{};]\s*$", script[max(0, i - 40):i]):
            # A regex literal: skip it so quotes inside it are not read as strings.
            j = i + 1
            in_class = False
            while j < n and (script[j] != "/" or in_class):
                if script[j] == "\\":
                    j += 1
                elif script[j] == "[":
                    in_class = True
                elif script[j] == "]":
                    in_class = False
                elif script[j] == "\n":
                    break
                j += 1
            i = j + 1
        else:
            i += 1


def decode_js(body: str) -> str:
    return body.replace("\\'", "'").replace('\\"', '"').replace("\\`", "`")


def encode_js(text: str, quote: str) -> str:
    text = text.replace("\\", "\\\\").replace("\r", "").replace("\n", " ")
    return text.replace(quote, "\\" + quote)


def _is_sentence(body: str) -> bool:
    escapes = set(re.findall(r"\\.", body))
    if not escapes <= ALLOWED_JS_ESCAPES:
        return False
    words = WORD_RE.findall(EXPRESSION_RE.sub(" ", re.sub(r"<[^>]+>", " ", body)))
    code = JS_CODE_HINTS.search(EXPRESSION_RE.sub("", body).strip())
    return len(words) >= MIN_WORDS_JS and " " in body and not code


# --- the list of text pieces ---------------------------------------------------

def extract_items(template: str) -> list[dict]:
    """Every piece of visible text in the template, in page order, with a stable id."""
    pieces, scripts = _html_pieces(template)
    raw = []
    section = "(page top)"
    for start, end, tag, page in pieces:
        span = _trim(template, start, end)
        if span is None or tag in MATH_TAGS:
            continue
        fragment = template[span[0]:span[1]]
        if tag in HEADING_TAGS:
            section = plain(fragment).strip()
        raw.append({"kind": "html", "start": span[0], "end": span[1], "page": page,
                    "section": section, "tag": tag, "text": fragment})
    for script_start, script_end in scripts:
        script = template[script_start:script_end]
        for start, end, quote in _js_strings(script, script_start):
            body = template[start + 1:end - 1]
            if not _is_sentence(body):
                continue
            before = list(JS_SECTION_RE.finditer(template, script_start, start))
            label = next((g for g in before[-1].groups() if g), "") if before else ""
            section = f"Tour card: {label}" if before and before[-1].group(1) else f"Script: {label}"
            raw.append({"kind": "js", "start": start, "end": end, "quote": quote, "page": "scripts",
                        "section": section, "tag": "script", "text": decode_js(body)})
    seen: dict[str, int] = {}
    for item in raw:
        key = item["kind"] + "\0" + item["text"]
        seen[key] = seen.get(key, 0) + 1
        item["id"] = hashlib.sha1(f"{key}\0{seen[key]}".encode()).hexdigest()[:10]
    return raw


def snapshot(items: list[dict]) -> list[dict]:
    keys = ("id", "kind", "page", "section", "tag", "text")
    return [{k: item[k] for k in keys} for item in items]


def write_snapshot(template: str) -> None:
    """Called by the generator: content/site_text.json = the text the site now shows."""
    SNAPSHOT_PATH.parent.mkdir(exist_ok=True)
    data = json.dumps(snapshot(extract_items(template)), ensure_ascii=False, indent=1)
    SNAPSHOT_PATH.write_text(data + "\n", encoding="utf-8")


# --- writing edits back ---------------------------------------------------------

FORMATTING_TAGS = {"b", "strong", "i", "em", "br"}  # Will may add or remove these freely


def _structure(fragment: str) -> list[str]:
    """The page elements inside a piece of text, other than bold/italic/line breaks."""
    tags = re.findall(r"<([a-zA-Z][\w-]*)([^>]*)>", fragment)
    # Name, class and link target: what the browser keeps when it re-saves edited HTML.
    keep = re.compile(r"""\b(class|href)\s*=\s*["']([^"']*)["']""")
    return sorted(name.lower() + str(sorted(keep.findall(attrs)))
                  for name, attrs in tags if name.lower() not in FORMATTING_TAGS)


def _problem(item: dict, new: str) -> str | None:
    """Why an edit can't go in safely, or None."""
    if sorted(PLACEHOLDER_RE.findall(item["text"])) != sorted(PLACEHOLDER_RE.findall(new)):
        return "a __NAME__ placeholder was changed or removed"
    if item["kind"] == "js" and EXPRESSION_RE.findall(item["text"]) != EXPRESSION_RE.findall(new):
        return "a ${...} part was changed or removed"
    if item["kind"] == "html" and _structure(item["text"]) != _structure(new):
        return "a link, label or other page element inside the text was removed or changed"
    prefix = CODE_PREFIX_RE.match(item["text"])
    if item["kind"] == "js" and prefix and not new.startswith(prefix.group(0)):
        return f"the code name {prefix.group(0)!r} at the start must stay"
    if '"""' in new:
        return 'it contains """'
    if not new.strip():
        return "the new text is empty"
    return None


def apply_edits(edits: list[dict]) -> tuple[list[str], list[str]]:
    """Write edits ({id, base, text}) into the generator. Returns (applied ids, problems)."""
    source, start = read_template()
    by_id = {item["id"]: item for item in extract_items(template_text(source, start))}
    applied, problems, changes = [], [], []
    for edit in edits:
        item = by_id.get(edit["id"])
        if item is None or item["text"] != edit["base"]:
            problems.append(f"{edit['id']}: the site's text here has changed since you edited it "
                            f"(was: {edit['base'][:60]!r})")
            continue
        reason = _problem(item, edit["text"])
        if reason:
            problems.append(f"{edit['id']}: {reason}")
            continue
        new = edit["text"].strip()
        if item["kind"] == "js":
            new = item["quote"] + encode_js(new, item["quote"]) + item["quote"]
        changes.append((item["start"], item["end"], new))
        applied.append(edit["id"])
    for item_start, item_end, new in sorted(changes, reverse=True):
        source = source[:start + item_start] + new + source[start + item_end:]
    GENERATOR_PATH.write_bytes(source.encode("utf-8"))
    return applied, problems


def check_scripts() -> bool:
    """Node syntax check of every inline script in index.html."""
    node = shutil.which("node")
    if node is None:
        print("Node.js not found: script syntax not checked")
        return False
    page = INDEX_PATH.read_text(encoding="utf-8")
    scripts = re.findall(r"<script(?![^>]*\bsrc=)(?![^>]*application/json)[^>]*>(.*?)</script>", page, re.S)
    ok = True
    for number, script in enumerate(scripts, 1):
        run = subprocess.run([node, "--check", "-"], input=script, capture_output=True, text=True,
                             encoding="utf-8")
        if run.returncode:
            ok = False
            print(f"Script {number} has a syntax error:\n{run.stderr[-800:]}")
    print(f"Checked {len(scripts)} scripts: {'OK' if ok else 'FAILED'}")
    return ok


def snapshot_hash() -> str:
    return hashlib.sha1(SNAPSHOT_PATH.read_bytes()).hexdigest()[:12]


def main(args: list[str]) -> int:
    command = args[0] if args else ""
    if command == "apply" and len(args) == 2:
        edits = json.loads(Path(args[1]).read_text(encoding="utf-8"))
        applied, problems = apply_edits(edits)
        print(json.dumps({"applied": applied, "problems": problems}, ensure_ascii=False, indent=1))
        return 0
    if command == "check":
        return 0 if check_scripts() else 1
    if command == "status":
        published = PUBLISHED_PATH.read_text().strip() if PUBLISHED_PATH.exists() else ""
        if published != snapshot_hash():
            print("The site text editor artifact is out of date: run the site-text skill "
                  "(pull Will's edits first, then republish the editor).", file=sys.stderr)
            return 2
        return 0
    if command == "hook":  # stdin: the Claude Code PostToolUse event for a Bash call
        event = json.loads(sys.stdin.read() or "{}")
        if "git commit" not in event.get("tool_input", {}).get("command", ""):
            return 0
        return main(["status"])
    if command == "mark-published":
        PUBLISHED_PATH.write_text(snapshot_hash() + "\n")
        return 0
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
