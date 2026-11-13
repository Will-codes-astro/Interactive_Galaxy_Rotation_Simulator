"""Read-only access to the SPARC mass-model table (Lelli, McGaugh & Schombert 2016).

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


def galaxy_rows(name: str, path: Path = SPARC_PATH) -> list[SparcRow]:
    """One galaxy's rows sorted by radius (name match ignores case)."""
    rows = sorted(
        (row for row in load_table(path) if row.name.upper() == name.upper()),
        key=lambda row: row.rad_kpc,
    )
    if len(rows) < 3:
        raise ValueError(f"No usable SPARC rows for {name} in {path}")
    return rows
