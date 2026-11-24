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
    python site_text.py layout OUT [BEFORE] text cut off / off screen at desktop and phone
                                           width; with BEFORE, prints only new problems
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
# Script strings like 'curvePlot:Watch the cyan line' start with a code name the page reads:
# the editor shows only the words after it, and apply puts the code name back.
CODE_PREFIX_RE = re.compile(r"^[a-z]\w*:(?! )")
MATH_TAGS = {"mi", "mo", "mn"}  # single maths symbols, not wording
ALLOWED_JS_ESCAPES = {"\\'", '\\"', "\\`"}

ATTR_RE = re.compile(r"""\s(title|data-tip|data-label|placeholder|alt|aria-label|content)\s*=\s*(["'])(.*?)\2""", re.S)
ATTR_LABELS = {"title": "Tooltip (shows on hover)", "data-tip": "Tooltip (shows on hover)",
               "data-label": "Label", "placeholder": "Hint inside an empty box",
               "alt": "Image description", "content": "Search-engine description",
               "aria-label": "Screen-reader label (read aloud, not shown)"}
TAG_LABELS = {"h1": "Heading", "h2": "Heading", "h3": "Heading", "h4": "Heading", "h5": "Small heading",
              "h6": "Small heading", "p": "Paragraph", "li": "List item", "button": "Button",
              "a": "Link", "label": "Control label", "dt": "Term", "dd": "Description",
              "td": "Table cell", "th": "Table heading", "figcaption": "Figure caption",
              "caption": "Table caption", "summary": "Click-to-open heading", "title": "Browser tab title",
              "option": "Menu choice", "legend": "Group heading", "blockquote": "Quote"}
TOUR_KEY_LABELS = {"title": "Card title", "text": "Card explanation", "action": "What to do",
                   "doneText": "Shown when the step is done", "tag": "Hint by the highlighted control",
                   "question": "Prediction question", "options": "Answer choice", "watch": "Watch hint",
                   "also": "Hint by another highlighted control", "label": "Tour name"}
PAGE_TITLES = {"head": "Browser tab and search results", "top bar": "Top bar and menu"}
TOUR_CARD_RE = re.compile(r"\{id:'([\w-]+)'(?:(?!\{id:).)*?\btitle:'((?:[^'\\\n]|\\.)*)'", re.S)
TOUR_NAME_RE = re.compile(r"(\w+):\{label:'((?:[^'\\\n]|\\.)*)'")
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
    """Return (pieces, tags, scripts): each run of inline text between block tags as
    (start, end, tag, page); every start tag as (start, end, name, page); and the
    (start, end) of every inline script."""
    scripts: list[tuple[int, int]] = []
    pieces, tags = [], []
    run_start = None
    block_tag = "body"
    div_depth = 0
    page = "head"
    page_depth = None
    pos = 0
    for match in TAG_RE.finditer(template):
        if match.start() < pos:
            continue
        closing, name, _ = match.group(1), (match.group(2) or "").lower(), match.group(3)
        if match.group(0).startswith("<!--"):
            continue
        if name == "body" and not closing:
            page = "top bar"
        if name in INLINE_TAGS:
            if not closing:
                tags.append((match.start(), match.end(), name, page))
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
                    page, page_depth = "shared", None
            else:
                div_depth += 1
                found = PAGE_RE.search(match.group(0))
                if found and 'class="page' in match.group(0):
                    page, page_depth = found.group(1), div_depth
        if not closing:
            block_tag = name
            tags.append((match.start(), match.end(), name, page))
        # Text right after this tag starts a new run.
        run_start = match.end()
    return pieces, tags, scripts


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


def _is_sentence(body: str, min_words: int = MIN_WORDS_JS) -> bool:
    escapes = set(re.findall(r"\\.", body))
    if not escapes <= ALLOWED_JS_ESCAPES:
        return False
    words = WORD_RE.findall(EXPRESSION_RE.sub(" ", re.sub(r"<[^>]+>", " ", body)))
    code = JS_CODE_HINTS.search(EXPRESSION_RE.sub("", body).strip())
    return len(words) >= min_words and (" " in body or min_words == 1) and not code


# --- the list of text pieces ---------------------------------------------------

def _humanise(name: str) -> str:
    words = re.sub(r"([a-z])([A-Z])", r"\1 \2", name).replace("_", " ").split()
    words = [w for w in words if w.lower() not in ("btn", "el")]
    return " ".join(words).capitalize()


def _attribute_items(template: str, tags) -> list[dict]:
    items = []
    for tag_start, tag_end, name, page in tags:
        tag = template[tag_start:tag_end]
        if name == "meta" and not re.search(r'(name|property)="[^"]*description"', tag):
            continue
        found = [(m.group(1), m.start(3) + tag_start, m.end(3) + tag_start, m.group(2), m.group(3))
                 for m in ATTR_RE.finditer(tag)]
        found = [f for f in found if len(re.findall(r"[A-Za-z]", f[4])) >= MIN_LETTERS_HTML]
        shown = [f for f in found if f[0] != "aria-label"]
        own_id = re.search(r'\sid="([^"]+)"', tag)
        after = plain(template[tag_end:template.find("<", tag_end)]).strip()
        before = plain(template[template.rfind(">", 0, tag_start) + 1:tag_start]).strip()
        where = (_humanise(own_id.group(1)) if own_id else
                 after if len(after) >= 3 else before if len(before) >= 3 else "")
        for attr, value_start, value_end, quote, value in found:
            linked = [(f[1], f[2]) for f in shown if f[0] != attr and f[4].strip() == value.strip()]
            if attr == "aria-label" and linked:
                continue  # it repeats a tooltip, and is updated with it
            if attr != "aria-label":
                linked += [(f[1], f[2]) for f in found
                           if f[0] == "aria-label" and f[4].strip() == value.strip()]
            items.append({"kind": "attr", "start": value_start, "end": value_end, "quote": quote,
                          "links": linked, "page": page, "label": ATTR_LABELS[attr],
                          "where": where, "text": html.unescape(value)})
    return items


def _block_end(script: str, open_at: int, strings: list[tuple[int, int, str]]) -> int:
    """Index just past the } that closes the { at open_at, skipping string literals."""
    depth, i, k = 0, open_at, 0
    while i < len(script):
        while k < len(strings) and strings[k][1] <= i:
            k += 1
        if k < len(strings) and strings[k][0] <= i:
            i = strings[k][1]
            continue
        depth += {"{": 1, "}": -1}.get(script[i], 0)
        i += 1
        if depth == 0:
            return i
    return len(script)


def _js_items(template: str, script_start: int, script_end: int) -> list[dict]:
    script = template[script_start:script_end]
    strings = [(a - script_start, b - script_start, q) for a, b, q in _js_strings(script, script_start)]
    tours = script.find("const TOURS={")
    tours_end = _block_end(script, tours + len("const TOURS="), strings) if tours >= 0 else -1
    card_titles = {m.group(1): decode_js(m.group(2)) for m in TOUR_CARD_RE.finditer(script)}
    tour_names = {m.group(1): decode_js(m.group(2)) for m in TOUR_NAME_RE.finditer(script)}
    items = []
    for a, b, quote in strings:
        body = script[a + 1:b - 1]
        key = re.search(r"(\w+)\s*:\s*\[?\s*$", script[max(0, a - 80):a])
        if not key and re.search(r",\s*$", script[max(0, a - 3):a]):
            key = re.search(r"(\w+):\[[^\]\n]*$", script[max(0, a - 300):a])
        key = key.group(1) if key else ""
        # Tour card fields are always wording, however short ('Inner star').
        tour_words = tours <= a < tours_end and key in TOUR_KEY_LABELS and key != "label"
        if not (_is_sentence(body, min_words=1 if tour_words else MIN_WORDS_JS)):
            continue
        text = decode_js(body)
        prefix = CODE_PREFIX_RE.match(text)
        prefix = prefix.group(0) if prefix else ""
        item = {"kind": "js", "start": script_start + a, "end": script_start + b, "quote": quote,
                "prefix": prefix, "text": text[len(prefix):]}
        if tours <= a < tours_end:
            card = list(re.finditer(r"\{id:'([\w-]+)'", script[tours:a]))
            tour = list(re.finditer(r"(\w+):\{label:", script[tours:a]))
            tour_name = tour_names.get(tour[-1].group(1), "") if tour else ""
            card_name = card_titles.get(card[-1].group(1), card[-1].group(1)) if card else "Tour names"
            item.update(page="tour", section=f"{tour_name}: {card_name}" if card else card_name,
                        label=TOUR_KEY_LABELS.get(key, "Tour message"))
        else:
            names = [g for m in JS_SECTION_RE.finditer(script, 0, a) for g in m.groups() if g]
            item.update(page="messages", section=_humanise(names[-1]) if names else "Page start",
                        label=TOUR_KEY_LABELS.get(key, "Message"))
        items.append(item)
    return items


def extract_items(template: str) -> list[dict]:
    """Every piece of text in the template (page text, tooltips and labels, script
    sentences), in page order, each with a stable id and a plain-words label."""
    pieces, tags, scripts = _html_pieces(template)
    page_items = []
    for start, end, tag, page in pieces:
        span = _trim(template, start, end)
        if span is None or tag in MATH_TAGS:
            continue
        page_items.append({"kind": "html", "start": span[0], "end": span[1], "page": page,
                           "tag": tag, "label": TAG_LABELS.get(tag, "Text"),
                           "text": template[span[0]:span[1]]})
    page_items += _attribute_items(template, tags)
    page_items.sort(key=lambda item: item["start"])
    section, page = "Top of the page", None
    for item in page_items:
        if item["page"] != page:  # a heading never carries over into the next page
            section, page = "Top of the page", item["page"]
        if item.get("tag") in HEADING_TAGS:
            section = plain(item["text"]).strip()
        item["section"] = PAGE_TITLES.get(item["page"], section)
    script_items = [i for a, b in scripts for i in _js_items(template, a, b)]
    raw = page_items + script_items
    seen: dict[str, int] = {}
    for item in raw:
        key = item["kind"] + "\0" + item["label"] + "\0" + item["text"]
        seen[key] = seen.get(key, 0) + 1
        item["id"] = hashlib.sha1(f"{key}\0{seen[key]}".encode()).hexdigest()[:10]
    return raw


def snapshot(items: list[dict]) -> list[dict]:
    keys = ("id", "kind", "page", "section", "label", "where", "text")
    return [{k: item[k] for k in keys if item.get(k)} for item in items]


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
    if item["kind"] == "html" and _structure(item["text"]) != _structure(new):
        return "a link, label or other page element inside the text was removed or changed"
    if item["kind"] != "html" and EXPRESSION_RE.findall(item["text"]) != EXPRESSION_RE.findall(new):
        return "a ${...} part was changed or removed"
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
        new = edit["text"].strip() if item["kind"] == "html" else re.sub(r"\s+", " ", edit["text"]).strip()
        if item["kind"] == "js":
            new = item["quote"] + encode_js(item["prefix"] + new, item["quote"]) + item["quote"]
        elif item["kind"] == "attr":
            new = new.replace("&", "&amp;").replace("<", "&lt;").replace(item["quote"], html.escape(item["quote"]))
            changes += [(a, b, new) for a, b in item["links"]]
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


LAYOUT_SIZES = ((1440, 900), (390, 844))  # desktop and phone, as in the build-log screenshots
LAYOUT_WAIT_MS = 15000  # virtual time for the page and its libraries to load
# Runs inside index.html: shows each page in turn and lists text that is cut off by its box
# or runs off the screen, then posts the list to the frame page around it.
LAYOUT_PROBE = r"""<script>
window.addEventListener('load',()=>setTimeout(()=>{
  const out=[],vw=document.documentElement.clientWidth;
  const scrolls=el=>{for(let a=el.parentElement;a;a=a.parentElement){const o=getComputedStyle(a).overflowX;if(o==='auto'||o==='scroll')return true;}return false;};
  const clipBox=el=>{for(let a=el.parentElement;a&&a!==document.body;a=a.parentElement){const c=getComputedStyle(a);if(/hidden|clip/.test(c.overflowX+c.overflowY))return a;}return null;};
  const check=page=>{
    const root=page?document.querySelector('.page[data-page="'+page+'"]'):document.querySelector('header,.topbar');
    if(!root)return;
    for(const el of root.querySelectorAll('*')){
      if(![...el.childNodes].some(n=>n.nodeType===3&&n.textContent.trim()))continue;
      const r=el.getBoundingClientRect(),cs=getComputedStyle(el);
      if(r.width<=1||r.height<=1||cs.visibility==='hidden'||el.closest('[hidden],[aria-hidden="true"]'))continue;
      let issue=null;
      if(/hidden|clip/.test(cs.overflowX)&&el.scrollWidth>el.clientWidth+1)issue='text cut off';
      else if(/hidden|clip/.test(cs.overflowY)&&el.scrollHeight>el.clientHeight+1)issue='text cut off at the bottom';
      else if(r.right>vw+1&&!scrolls(el))issue='runs off the screen';
      else{const b=clipBox(el);if(b){const q=b.getBoundingClientRect();if(q.width>1&&q.height>1&&((r.left<q.right&&r.right>q.right+1)||(r.top<q.bottom&&r.bottom>q.bottom+1)))issue='spills out of its box';}}
      if(issue)out.push([page||'top bar',issue,el.tagName.toLowerCase()+(el.id?'#'+el.id:'')+(el.classList[0]?'.'+el.classList[0]:''),el.textContent.trim().replace(/\s+/g,' ').slice(0,70)].join(' | '));
    }
  };
  check(null);
  for(const p of [...document.querySelectorAll('.page')].map(e=>e.dataset.page)){if(typeof showPage==='function')showPage(p,{scroll:false});check(p);}
  parent.postMessage('layout:'+JSON.stringify(out),'*');
},1500));
</script>"""


def layout_issues() -> list[str]:
    """Text that is cut off or runs off the screen on any page of index.html, at each size."""
    browser = next((b for b in (shutil.which("chrome"), shutil.which("msedge"),
                                r"C:\Program Files\Google\Chrome\Application\chrome.exe",
                                r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")
                    if b and Path(b).exists()), None)
    if browser is None:
        raise SystemExit("Chrome or Edge not found: layout not checked")
    page = INDEX_PATH.read_text(encoding="utf-8")
    probe = BASE_DIR / "_layout_probe.html"  # next to index.html so relative files still load
    probe.write_text(page.replace("</body>", LAYOUT_PROBE + "</body>"), encoding="utf-8")
    # Headless Chrome won't make a window narrower than ~500 px, so the page runs in a frame
    # of the exact width; the frame page writes the probe's report where --dump-dom reads it.
    frame = BASE_DIR / "_layout_frame.html"
    issues = []
    try:
        for width, height in LAYOUT_SIZES:
            frame.write_text(
                f'<body style="margin:0"><iframe src="{probe.name}" style="width:{width}px;height:{height}px;'
                'border:0"></iframe><script>addEventListener("message",e=>{if(String(e.data).startsWith("layout:"))'
                '{const p=document.createElement("pre");p.id="layoutReport";p.textContent=e.data.slice(7);'
                'document.body.append(p)}})</script></body>', encoding="utf-8")
            run = subprocess.run([browser, "--headless=new", "--use-angle=swiftshader", "--enable-unsafe-swiftshader",
                                  "--hide-scrollbars", f"--window-size={width + 40},{height + 40}",
                                  f"--virtual-time-budget={LAYOUT_WAIT_MS}", "--dump-dom", frame.as_uri()],
                                 capture_output=True, text=True, encoding="utf-8", timeout=120)
            found = re.search(r'<pre id="layoutReport">(.*?)</pre>', run.stdout, re.S)
            if not found:
                raise SystemExit(f"The layout probe did not report at {width}px:\n{run.stderr[-600:]}")
            issues += [f"{width}px | {line}" for line in json.loads(html.unescape(found.group(1)))]
    finally:
        probe.unlink(missing_ok=True)
        frame.unlink(missing_ok=True)
    return sorted(set(issues))


def snapshot_hash() -> str:
    return hashlib.sha1(SNAPSHOT_PATH.read_bytes()).hexdigest()[:12]


def main(args: list[str]) -> int:
    command = args[0] if args else ""
    if command == "apply" and len(args) == 2:
        edits = json.loads(Path(args[1]).read_text(encoding="utf-8"))
        applied, problems = apply_edits(edits)
        print(json.dumps({"applied": applied, "problems": problems}, ensure_ascii=False, indent=1))
        return 0
    if command == "layout" and len(args) in (2, 3):
        issues = layout_issues()
        Path(args[1]).write_text(json.dumps(issues, ensure_ascii=False, indent=1), encoding="utf-8")
        if len(args) == 3:
            before = set(json.loads(Path(args[2]).read_text(encoding="utf-8")))
            new = [i for i in issues if i not in before]
            print(chr(10).join(new) if new else "No new layout problems")
            return 1 if new else 0
        print(f"{len(issues)} layout notes saved")
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
