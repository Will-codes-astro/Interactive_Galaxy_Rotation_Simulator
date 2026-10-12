"""Generate high-resolution PNG and SVG QR codes for the EPQ dashboard."""

from __future__ import annotations

from pathlib import Path

import qrcode
from qrcode.constants import ERROR_CORRECT_H
from qrcode.image.svg import SvgPathImage


BASE_DIR = Path(__file__).resolve().parent
OUTPUT_PNG = BASE_DIR / "epq_qr_code.png"
OUTPUT_SVG = BASE_DIR / "epq_qr_code.svg"
TARGET_URL = "https://will-codes-astro.github.io/Interactive_Galaxy_Rotation_Simulator/"


def main() -> None:
    code = qrcode.QRCode(
        version=None,
        error_correction=ERROR_CORRECT_H,
        box_size=24,
        border=5,
    )
    code.add_data(TARGET_URL)
    code.make(fit=True)

    png = code.make_image(fill_color="#07111d", back_color="#ffffff")
    png.save(OUTPUT_PNG)
    svg = code.make_image(
        image_factory=SvgPathImage,
        fill_color="#07111d",
        back_color="#ffffff",
    )
    svg.save(OUTPUT_SVG)

    print(f"QR target: {TARGET_URL}")
    print(f"Generated {OUTPUT_PNG.name} ({OUTPUT_PNG.stat().st_size:,} bytes)")
    print(f"Generated {OUTPUT_SVG.name} ({OUTPUT_SVG.stat().st_size:,} bytes)")
    print("QR code ready for the configured GitHub Pages deployment.")


if __name__ == "__main__":
    main()