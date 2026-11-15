"""Read-only access to the SPARC tables (Lelli, McGaugh & Schombert 2016): the mass models
(table2.dat) and the galaxy sample (table1.mrt, for distances and inclinations).

Every script that needs SPARC data goes through this module, so the raw file is
checked the same way everywhere. It never writes, downloads or replaces the file.
Column meanings, units and provenance: docs/DATA_DICTIONARY.md.

Standard library only, so the site generator can use it without extra packages.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import NamedTuple

ROOT = Path(__file__).resolve().parent.parent
SPARC_PATH = ROOT / "data" / "sparc" / "table2.dat"

# SHA-256 of table2.dat with line endings normalised to LF (git/OneDrive may store CRLF).
# Recorded 2026-11-13 from the copy downloaded from VizieR on 2026-09-25.
SPARC_SHA256 = "40144d85a2278064c281bbca7fca2a614098b43002e7d4f43f32ecd7f9728995"
HEADER_NAMES = ("Name", "Dist", "Rad", "Vobs", "e_Vobs", "Vgas", "Vdisk", "Vbulge", "SBdisk", "SBbulge")

# Allowed ranges from the VizieR column descriptions in the file header.
VALID_RANGES = {
    "dist_mpc": (0.98, 127.8),
    "rad_kpc": (0.08, 108.31),
    "v_obs_kms": (1.4, 383.0),
    "e_v_obs_kms": (0.2, 75.0),
    "v_gas_kms": (-17.0, 87.0),
    "v_disk_kms": (1.0, 372.0),
    "v_bulge_kms": (0.0, 390.0),
}
EXPECTED_GALAXIES = 175

# SPARC Table 1 (galaxy sample): distances and inclinations with their errors.
TABLE1_PATH = ROOT / "data" / "sparc" / "table1.mrt"
# SHA-256 with line endings normalised to LF. Recorded 2026-11-15 from the copy downloaded that day
# from astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt (server date 28 Mar 2016).
TABLE1_SHA256 = "5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243"
# Positions of the fields used, counting whitespace-separated fields (Galaxy T D e_D f_D Inc e_Inc ...).
# The data lines are offset from the file's own byte-by-byte description, so they are split on spaces;
# no field used here contains a space.
TABLE1_FIELDS = {"name": 0, "dist_mpc": 2, "e_dist_mpc": 3, "inc_deg": 5, "e_inc_deg": 6, "hubble_type": 1,
                 "lum_36_1e9_lsun": 7, "hi_mass_1e9_msun": 13, "v_flat_kms": 15}
TABLE1_FIELD_COUNT = 19


class SparcRow(NamedTuple):
    """One radius of one galaxy, exactly as published (no M/L ratios applied)."""

    name: str
    dist_mpc: float  # assumed distance, Mpc
    rad_kpc: float  # galactocentric radius, kpc
    v_obs_kms: float  # observed circular speed, km/s
    e_v_obs_kms: float  # 1-sigma uncertainty of v_obs, km/s
    v_gas_kms: float  # gas contribution, km/s; can be negative (see DATA_DICTIONARY.md)
    v_disk_kms: float  # stellar disk contribution for M/L = 1, km/s
    v_bulge_kms: float  # bulge contribution for M/L = 1, km/s
    sb_disk: float  # disk surface brightness, Lsun/pc^2
    sb_bulge: float  # bulge surface brightness, Lsun/pc^2


def file_sha256(path: Path = SPARC_PATH) -> str:
    """SHA-256 of the file with CRLF normalised to LF."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_table(path: Path = SPARC_PATH, check_hash: bool = True) -> list[SparcRow]:
    """Parse and validate the whole table. Raises ValueError on any problem."""
    if not path.is_file():
        raise FileNotFoundError(f"SPARC table missing: {path}. See docs/DATA_DICTIONARY.md.")
    if check_hash and path.resolve() == SPARC_PATH.resolve():
        digest = file_sha256(path)
        if digest != SPARC_SHA256:
            raise ValueError(
                f"{path} has changed (sha256 {digest}). Raw SPARC data must not be edited; "
                "restore it with git."
            )

    lines = path.read_text(encoding="utf-8").splitlines()
    # VizieR TSV layout: '#' comment lines, then a column-name line, a unit line and a dashed line.
    data_lines = [(n, line) for n, line in enumerate(lines, 1) if line.strip() and not line.startswith("#")]
    if len(data_lines) < 3 or data_lines[0][1].split() != list(HEADER_NAMES):
        raise ValueError(f"{path}: column header is not {' '.join(HEADER_NAMES)}")
    rows = []
    for line_number, line in data_lines[3:]:
        parts = line.split()
        if len(parts) != 10:
            raise ValueError(f"{path}:{line_number}: expected 10 columns, found {len(parts)}")
        try:
            row = SparcRow(parts[0], *map(float, parts[1:]))
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: {error}") from None
        for field, (low, high) in VALID_RANGES.items():
            value = getattr(row, field)
            if not low <= value <= high:
                raise ValueError(f"{path}:{line_number}: {field}={value} outside [{low}, {high}]")
        rows.append(row)

    galaxies = {row.name for row in rows}
    if len(galaxies) != EXPECTED_GALAXIES:
        raise ValueError(f"{path}: expected {EXPECTED_GALAXIES} galaxies, found {len(galaxies)}")
    return rows


class GalaxyInfo(NamedTuple):
    """One galaxy's distance and inclination with their published errors, and its catalogue properties (SPARC Table 1)."""

    name: str
    dist_mpc: float  # assumed distance, Mpc (the same value as table2's Dist)
    e_dist_mpc: float  # mean error on the distance, Mpc
    inc_deg: float  # disk inclination, degrees (90 = edge-on)
    e_inc_deg: float  # mean error on the inclination, degrees
    hubble_type: float  # T: 0 = S0, 1 = Sa ... 5 = Sc ... 10 = Im, 11 = BCD (Table 1 note 1; a whole number)
    lum_36_1e9_lsun: float  # L[3.6]: total luminosity at 3.6 um, 10^9 Lsun
    hi_mass_1e9_msun: float  # MHI: total atomic hydrogen mass, 10^9 Msun
    v_flat_kms: float  # Vflat: asymptotically flat rotation speed, km/s (0 where SPARC gives none)


def load_galaxy_info(path: Path = TABLE1_PATH) -> dict[str, GalaxyInfo]:
    """Parse SPARC Table 1 (machine-readable text), keyed by galaxy name. Raises ValueError on any problem."""
    if not path.is_file():
        raise FileNotFoundError(f"SPARC Table 1 missing: {path}. See docs/DATA_DICTIONARY.md.")
    if path.resolve() == TABLE1_PATH.resolve() and file_sha256(path) != TABLE1_SHA256:
        raise ValueError(f"{path} has changed. Raw SPARC data must not be edited; restore it with git.")
    lines = path.read_text(encoding="ascii").splitlines()
    # The byte-by-byte description ends with the last dashed line; one galaxy per line follows.
    last_rule = max(n for n, line in enumerate(lines) if line.startswith("-" * 20))
    info = {}
    for line_number, line in enumerate(lines[last_rule + 1:], last_rule + 2):
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) != TABLE1_FIELD_COUNT:
            raise ValueError(f"{path}:{line_number}: expected {TABLE1_FIELD_COUNT} fields, found {len(parts)}")
        try:
            name = parts[TABLE1_FIELDS["name"]]
            values = [float(parts[TABLE1_FIELDS[field]]) for field in GalaxyInfo._fields[1:]]
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: {error}") from None
        row = GalaxyInfo(name, *values)
        if not (row.dist_mpc > 0 and row.e_dist_mpc >= 0 and 0 < row.inc_deg <= 90 and row.e_inc_deg >= 0):
            raise ValueError(f"{path}:{line_number}: distance or inclination out of range: {row}")
        if not (row.hubble_type in range(12) and row.lum_36_1e9_lsun >= 0 and row.hi_mass_1e9_msun >= 0):
            raise ValueError(f"{path}:{line_number}: Hubble type, luminosity or HI mass out of range: {row}")
        info[name.upper()] = row
    if len(info) != EXPECTED_GALAXIES:
        raise ValueError(f"{path}: expected {EXPECTED_GALAXIES} galaxies, found {len(info)}")
    return info


def galaxy_rows(name: str, path: Path = SPARC_PATH) -> list[SparcRow]:
    """One galaxy's rows sorted by radius (name match ignores case)."""
    rows = sorted(
        (row for row in load_table(path) if row.name.upper() == name.upper()),
        key=lambda row: row.rad_kpc,
    )
    if len(rows) < 3:
        raise ValueError(f"No usable SPARC rows for {name} in {path}")
    return rows
