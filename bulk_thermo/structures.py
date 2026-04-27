"""Helpers to build or load periodic crystal structures.

This module isolates all logic related to constructing ASE ``Atoms`` objects
for bulk workflows. The rest of the code can then work with a ready-to-use
structure without caring whether it came from a geometry file or from lattice
parameters and fractional coordinates.
"""

from collections import Counter
from functools import reduce
from math import gcd

from ase.io import read
from ase.spacegroup import crystal


def setup_crystal(lattice_geometry, lattice_basis, calculator):
    """Build a periodic crystal from lattice parameters and fractional basis."""
    atom_symbols = [atom["symbol"] for atom in lattice_basis]
    atom_positions = [atom["position"] for atom in lattice_basis]

    atoms = crystal(
        atom_symbols,
        basis=atom_positions,
        spacegroup=lattice_geometry["Inum"],
        cellpar=[
            lattice_geometry["a"],
            lattice_geometry["b"],
            lattice_geometry["c"],
            lattice_geometry["alpha"],
            lattice_geometry["beta"],
            lattice_geometry["gamma"],
        ],
        pbc=(1, 1, 1),
    )

    atoms.calc = calculator
    return atoms


def setup_structure_from_file(geometry_file, calculator):
    """Read a periodic structure file and attach the calculator."""
    atoms = read(geometry_file)
    atoms.calc = calculator
    return atoms


def build_structure(args, calculator):
    """Build the structure either from file or from lattice JSON data."""
    if args.geometry_file is not None:
        return setup_structure_from_file(args.geometry_file, calculator)

    if args.lattice_geometry is None or args.lattice_basis is None:
        raise ValueError(
            "Provide either --geometry-file or both --lattice-geometry and --lattice-basis."
        )

    return setup_crystal(args.lattice_geometry, args.lattice_basis, calculator)


def infer_formula_units(atoms):
    """Infer the number of reduced formula units in the simulation cell.

    The reduced number of formula units is estimated as the greatest common
    divisor of the elemental populations in the cell.
    """
    counts = Counter(atoms.get_chemical_symbols())
    return reduce(gcd, counts.values())


def formula_strings(atoms, formula_units):
    """Return full-cell formula and reduced formula-unit string.

    The element order is preserved from the input structure so that the output
    remains familiar to the user.
    """
    symbols = atoms.get_chemical_symbols()
    counts = Counter(symbols)
    ordered_symbols = list(dict.fromkeys(symbols))

    def fmt(divisor):
        parts = []
        for symbol in ordered_symbols:
            value = counts[symbol] // divisor
            parts.append(symbol if value == 1 else f"{symbol}{value}")
        return "".join(parts)

    return fmt(1), fmt(formula_units)
