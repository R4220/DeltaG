"""Phonon utilities for bulk calculations.

The functions here wrap ASE's finite-displacement phonon machinery and keep
all phonon-specific file handling away from the higher-level thermodynamic
workflow.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from ase.phonons import Phonons


def run_phonons(atoms, calculator, supercell, delta, dos_kpts, bandpath=None, output_dir="."):
    """Compute finite-displacement phonons, DOS, and optional band structure."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ph = Phonons(
        atoms,
        calculator,
        supercell=tuple(supercell),
        delta=delta,
        name=str(output_dir / "phonon"),
    )

    try:
        ph.run()
        ph.read(acoustic=True)
    finally:
        # Remove displacement-force directories even if something fails later.
        ph.clean()

    if bandpath is not None:
        path = atoms.cell.bandpath(bandpath, npoints=100)
        band_structure = ph.get_band_structure(path)
    else:
        band_structure = None

    dos = ph.get_dos(kpts=tuple(dos_kpts)).sample_grid(
        npts=3000,
        width=5e-4,
        xmin=0.0,
    )

    return ph, dos, band_structure


def plot_phonon_bs_and_dos(dos, band_structure=None, output_path="phonon_BS_and_DOS.png", emax=0.035):
    """Plot phonon band structure and DOS in a compact figure."""
    fig = plt.figure(figsize=(7, 4))
    ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])

    if band_structure is not None:
        band_structure.plot(ax=ax, emin=0.0, emax=emax)
    else:
        ax.set_ylim(0.0, emax)
        ax.set_ylabel("Energy (eV)")
        ax.set_xticks([])
        ax.set_title("Phonon DOS only")

    dosax = fig.add_axes([0.8, 0.07, 0.17, 0.85])
    dosax.fill_between(
        dos.get_weights(),
        dos.get_energies(),
        y2=0,
        color="grey",
        edgecolor="k",
        lw=1,
    )

    dosax.set_ylim(0.0, emax)
    dosax.set_yticks([])
    dosax.set_xticks([])
    dosax.set_xlabel("DOS", fontsize=14)

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
