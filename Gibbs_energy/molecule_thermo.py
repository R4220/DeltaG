"""Ideal-gas molecular thermochemistry with ASE and MACE.

This module is dedicated to isolated molecules in the gas phase. The workflow
implemented here is:

1. read and box the molecule;
2. relax the geometry;
3. compute vibrational frequencies;
4. evaluate ideal-gas thermodynamic quantities with ASE.
"""

from collections import Counter
from pathlib import Path

from ase import units
from ase.io import read, write
from ase.optimize import QuasiNewton
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations

from .calculators import setup_calculator
from .io_utils import write_output


def setup_molecule(fname_geometry, vacuum, calculator):
    """Read a molecule, add a vacuum box, disable PBC, and attach calculator."""
    atoms = read(fname_geometry)

    dimensions = atoms.positions.max(axis=0) - atoms.positions.min(axis=0)
    atoms.cell = dimensions + vacuum
    atoms.center()
    atoms.pbc = False
    atoms.calc = calculator

    return atoms


def formula_string(atoms):
    """Return molecular formula preserving the original element order."""
    symbols = atoms.get_chemical_symbols()
    counts = Counter(symbols)
    ordered_symbols = list(dict.fromkeys(symbols))

    parts = []
    for symbol in ordered_symbols:
        value = counts[symbol]
        parts.append(symbol if value == 1 else f"{symbol}{value}")

    return "".join(parts)


def run_molecular_vibrations(atoms, output_dir):
    """Run finite-displacement vibrational analysis for an isolated molecule."""
    output_dir = Path(output_dir)
    vib_name = str(output_dir / "vib")

    vib = Vibrations(atoms, name=vib_name)
    vib.run()

    vib_energies = vib.get_energies()
    vib.summary(log=str(output_dir / "vibrations_summary.txt"))

    return vib, vib_energies


def compute_molecule_thermo(
    atoms,
    vib_energies,
    potentialenergy,
    mol_geometry,
    symmetry_number,
    spin,
    temperature,
    pressure,
):
    """Compute ideal-gas molecular U, H, S, F and G at one temperature and pressure."""
    thermo = IdealGasThermo(
        vib_energies=vib_energies,
        potentialenergy=potentialenergy,
        atoms=atoms,
        geometry=mol_geometry,
        symmetrynumber=symmetry_number,
        spin=spin,
    )

    H = thermo.get_enthalpy(temperature=temperature)
    U = H - units.kB * temperature
    S = thermo.get_entropy(
        temperature=temperature,
        pressure=pressure,
        verbose=False,
    )
    G = thermo.get_gibbs_energy(
        temperature=temperature,
        pressure=pressure,
    )

    TS = temperature * S
    F = U - TS

    return {
        "U": U,
        "H": H,
        "S": S,
        "TS": TS,
        "F": F,
        "G": G,
    }


def run_molecule(args):
    """Run full ideal-gas molecular thermochemistry workflow."""
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    calc = setup_calculator(args.model_path, device=args.device)

    # Prepare an isolated molecule in a finite vacuum box before relaxation.
    atoms = setup_molecule(
        fname_geometry=args.geometry_file,
        vacuum=args.vacuum,
        calculator=calc,
    )

    write(output_dir / "initial_molecule.xyz", atoms)

    optimizer = QuasiNewton(
        atoms,
        logfile=str(output_dir / "relax_molecule.log"),
        trajectory=str(output_dir / "relax_molecule.traj"),
    )
    optimizer.run(fmax=args.fmax)

    write(output_dir / "relaxed_molecule.xyz", atoms)

    potentialenergy = atoms.get_potential_energy()

    vib, vib_energies = run_molecular_vibrations(
        atoms=atoms,
        output_dir=output_dir,
    )

    # Thermodynamic quantities are evaluated at a single user-selected T and p.
    result = compute_molecule_thermo(
        atoms=atoms,
        vib_energies=vib_energies,
        potentialenergy=potentialenergy,
        mol_geometry=args.mol_geometry,
        symmetry_number=args.symmetry_number,
        spin=args.spin,
        temperature=args.temperature,
        pressure=args.pressure,
    )

    formula = formula_string(atoms)

    output_lines = [
        "# Ideal-gas molecular thermochemistry report",
        "",
        "[System]",
        f"Geometry file              : {args.geometry_file}",
        f"Molecular formula          : {formula}",
        f"Molecular geometry         : {args.mol_geometry}",
        f"Symmetry number            : {args.symmetry_number}",
        f"Total electronic spin      : {args.spin}",
        "",
        "[Calculation settings]",
        f"Model path                 : {args.model_path}",
        f"Device                     : {args.device}",
        f"Temperature (K)            : {args.temperature}",
        f"Pressure (Pa)              : {args.pressure}",
        f"Vacuum size (A)            : {args.vacuum}",
        f"Force threshold (eV/A)     : {args.fmax}",
        "",
        "[Meaning of the quantities]",
        "E_pot = static potential energy of the relaxed molecule",
        "U(T)  = internal energy including translational, rotational, and vibrational terms",
        "H(T)  = enthalpy in the ideal-gas approximation",
        "F(T)  = Helmholtz free energy = U - T*S",
        "G(T,p)= Gibbs free energy = H - T*S",
        "",
        "[Results per molecule]",
        f"E_pot (eV)                 : {potentialenergy:.16f}",
        f"U(T) (eV)                  : {result['U']:.16f}",
        f"S(T,p) (eV/K)              : {result['S']:.16f}",
        f"T*S (eV)                   : {result['TS']:.16f}",
        f"H(T) (eV)                  : {result['H']:.16f}",
        f"F(T) (eV)                  : {result['F']:.16f}",
        f"G(T,p) (eV)                : {result['G']:.16f}",
        "",
        "[Generated files]",
        "Initial molecule           : initial_molecule.xyz",
        "Relaxed molecule           : relaxed_molecule.xyz",
        "Relax log                  : relax_molecule.log",
        "Relax trajectory           : relax_molecule.traj",
        "Vibrational summary        : vibrations_summary.txt",
        "",
        "[Important note]",
        "This uses IdealGasThermo, so it is appropriate for isolated gas-phase molecules.",
        "For adsorbed molecules, do not use ideal-gas translational/rotational terms.",
    ]

    write_output(output_dir / "molecule_thermo_summary.out", output_lines)

    # Optionally remove the finite-displacement files after extracting results.
    if args.clean_vib:
        vib.clean()

    print(f"Molecular thermodynamic summary written to: {output_dir / 'molecule_thermo_summary.out'}")
