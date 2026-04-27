"""Lightweight helpers for human-readable output files.

The workflows in this directory generate summary reports and simple tabulated
data. These helpers keep the formatting logic in one place so the scientific
modules can focus on the actual calculations.
"""

from pathlib import Path


def write_output(output_path, lines):
    """Write a plain-text report, creating parent directories if needed."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_table(path, header, temperatures, values):
    """Write a simple temperature/value table with a custom header."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as handle:
        handle.write(header)
        for temperature, value in zip(temperatures, values):
            handle.write(f"{temperature:22.12f} {value:22.12f}\n")
