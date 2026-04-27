"""Harmonic thermochemistry for finite patches in vacuum.

This workflow targets finite fragments such as surface patches or supported
clusters represented in a large vacuum box. Unlike the molecular workflow, it
does not add ideal-gas translational or rotational terms. Instead, it uses the
harmonic approximation through ASE's ``HarmonicThermo`` class.
"""

from collections import Counter
from pathlib import Path

import numpy as np
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.optimize import QuasiNewton
from ase.thermochemistry import HarmonicThermo
from ase.vibrations import Vibrations

from .calculators import setup_calculator
from .io_utils import write_output


def setup_patch(geometry_file, vacuum, calculator):
    """Read a finite structure, box it in vacuum, disable PBC, and attach a calculator."""
    atoms = read(geometry_file)

    dimensions = atoms.positions.max(axis=0) - atoms.positions.min(axis=0)
    atoms.cell = dimensions + vacuum
    atoms.center()
    atoms.pbc = False
    atoms.calc = calculator

    return atoms


def formula_string(atoms):
    """Return a chemical formula preserving the original element order."""
    symbols = atoms.get_chemical_symbols()
    counts = Counter(symbols)
    ordered_symbols = list(dict.fromkeys(symbols))

    parts = []
    for symbol in ordered_symbols:
        value = counts[symbol]
        parts.append(symbol if value == 1 else f"{symbol}{value}")

    return "".join(parts)


def normalize_atom_indices(indices, natoms, label):
    """Validate a user-provided list of 0-based atom indices."""
    if indices is None:
        return None

    normalized = [int(index) for index in indices]
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{label} must not contain duplicates.")

    invalid = [index for index in normalized if index < 0 or index >= natoms]
    if invalid:
        invalid_list = ", ".join(str(index) for index in invalid)
        raise ValueError(
            f"{label} contains out-of-range atom indices for a structure with {natoms} atoms: {invalid_list}"
        )

    return normalized


def add_fixed_atom_constraints(atoms, fixed_indices):
    """Apply optional fixed-atom constraints without discarding pre-existing ones."""
    if not fixed_indices:
        return

    existing = atoms.constraints
    if not existing:
        constraints = []
    elif isinstance(existing, list):
        constraints = list(existing)
    elif isinstance(existing, tuple):
        constraints = list(existing)
    else:
        constraints = [existing]

    constraints.append(FixAtoms(indices=fixed_indices))
    atoms.set_constraint(constraints)


def run_patch_vibrations(atoms, output_dir, delta, vibration_indices):
    """Run a finite-displacement vibrational analysis on a finite patch."""
    output_dir = Path(output_dir)
    vib_name = str(output_dir / "vib")

    vib = Vibrations(
        atoms,
        indices=vibration_indices,
        name=vib_name,
        delta=delta,
    )
    vib.run()

    vib_energies = vib.get_energies()
    vib.summary(log=str(output_dir / "vibrations_summary.txt"))

    return vib, vib_energies


def compute_patch_thermo_at_temperature(thermo, potentialenergy, temperature):
    """Compute harmonic thermodynamic quantities for a finite patch."""
    e_zpe = thermo.get_ZPE_correction()

    if temperature <= 0.0:
        e_thermal_vib = 0.0
        U = potentialenergy + e_zpe
        S = 0.0
        TS = 0.0
        F = U
    else:
        U = thermo.get_internal_energy(temperature=temperature, verbose=False)
        S = thermo.get_entropy(temperature=temperature, verbose=False)
        F = thermo.get_helmholtz_energy(temperature=temperature, verbose=False)
        e_thermal_vib = U - potentialenergy - e_zpe
        TS = temperature * S

    return {
        "E_pot": potentialenergy,
        "E_ZPE": e_zpe,
        "E_thermal_vib": e_thermal_vib,
        "U": U,
        "S": S,
        "TS": TS,
        "F": F,
        "G_approx": F,
    }


def write_patch_temperature_table(output_path, temperatures, thermo, potentialenergy):
    """Write harmonic patch thermodynamics as a function of temperature."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# T[K] E_pot[eV] E_ZPE[eV] E_thermal_vib[eV] U[eV] "
            "S[eV/K] TS[eV] F[eV] G_approx[eV]\n"
        )

        for temperature in temperatures:
            result = compute_patch_thermo_at_temperature(
                thermo=thermo,
                potentialenergy=potentialenergy,
                temperature=temperature,
            )

            handle.write(
                f"{temperature:16.6f} "
                f"{result['E_pot']:18.10f} "
                f"{result['E_ZPE']:18.10f} "
                f"{result['E_thermal_vib']:18.10f} "
                f"{result['U']:18.10f} "
                f"{result['S']:18.10f} "
                f"{result['TS']:18.10f} "
                f"{result['F']:18.10f} "
                f"{result['G_approx']:18.10f}\n"
            )


def run_patch(args):
    """Run full harmonic thermochemistry for a finite patch in vacuum."""
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    calc = setup_calculator(args.model_path, device=args.device)
    atoms = setup_patch(
        geometry_file=args.geometry_file,
        vacuum=args.vacuum,
        calculator=calc,
    )

    natoms = len(atoms)
    fixed_indices = normalize_atom_indices(args.fixed_indices, natoms, "fixed_indices")
    vibration_indices = normalize_atom_indices(args.vibration_indices, natoms, "vibration_indices")

    if args.t_step <= 0.0:
        raise ValueError("t_step must be positive.")
    if args.t_max < args.t_min:
        raise ValueError("t_max must be greater than or equal to t_min.")

    if fixed_indices is not None and len(fixed_indices) == natoms:
        raise ValueError("fixed_indices cannot include every atom in the patch.")
    if vibration_indices == []:
        raise ValueError("vibration_indices cannot be empty.")

    if fixed_indices and vibration_indices is not None:
        overlap = sorted(set(fixed_indices) & set(vibration_indices))
        if overlap:
            overlap_list = ", ".join(str(index) for index in overlap)
            raise ValueError(
                f"vibration_indices must not overlap fixed_indices. Overlapping indices: {overlap_list}"
            )

    add_fixed_atom_constraints(atoms, fixed_indices)

    write(output_dir / "initial_patch.xyz", atoms)

    optimizer = QuasiNewton(
        atoms,
        logfile=str(output_dir / "relax_patch.log"),
        trajectory=str(output_dir / "relax_patch.traj"),
    )
    optimizer.run(fmax=args.fmax)

    write(output_dir / "relaxed_patch.xyz", atoms)

    potentialenergy = atoms.get_potential_energy()

    vib, vib_energies = run_patch_vibrations(
        atoms=atoms,
        output_dir=output_dir,
        delta=args.delta,
        vibration_indices=vibration_indices,
    )

    thermo = HarmonicThermo(
        vib_energies=vib_energies,
        potentialenergy=potentialenergy,
    )

    temperatures = np.arange(args.t_min, args.t_max + args.t_step, args.t_step)
    write_patch_temperature_table(
        output_path=output_dir / "patch_thermo_temperature.dat",
        temperatures=temperatures,
        thermo=thermo,
        potentialenergy=potentialenergy,
    )

    result = compute_patch_thermo_at_temperature(
        thermo=thermo,
        potentialenergy=potentialenergy,
        temperature=args.temperature,
    )

    formula = formula_string(atoms)
    fixed_label = tuple(fixed_indices) if fixed_indices else "None"
    if vibration_indices is None:
        vibration_label = "Automatic (all unconstrained atoms)"
    else:
        vibration_label = tuple(vibration_indices)

    output_lines = [
        "# Harmonic thermochemistry report for a finite patch",
        "",
        "[System]",
        f"Geometry file              : {args.geometry_file}",
        f"Patch formula              : {formula}",
        f"Number of atoms            : {natoms}",
        f"Fixed atom indices         : {fixed_label}",
        f"Vibrated atom indices      : {vibration_label}",
        "",
        "[Calculation settings]",
        f"Model path                 : {args.model_path}",
        f"Device                     : {args.device}",
        f"Temperature selected (K)   : {args.temperature}",
        f"Vacuum size (A)            : {args.vacuum}",
        f"Force threshold (eV/A)     : {args.fmax}",
        f"Vibrational delta (A)      : {args.delta}",
        "",
        "[Meaning of the quantities]",
        "E_pot = static potential energy of the relaxed patch",
        "U(T)  = harmonic internal energy including ZPE and thermal vibrations",
        "F(T)  = Helmholtz free energy in the harmonic approximation",
        "G(T)  = approximated here as F(T), assuming negligible pV",
        "",
        "[Results per patch]",
        f"E_pot (eV)                 : {result['E_pot']:.16f}",
        f"E_ZPE (eV)                 : {result['E_ZPE']:.16f}",
        f"E_thermal_vib (eV)         : {result['E_thermal_vib']:.16f}",
        f"U(T) (eV)                  : {result['U']:.16f}",
        f"S(T) (eV/K)                : {result['S']:.16f}",
        f"T*S (eV)                   : {result['TS']:.16f}",
        f"F(T) (eV)                  : {result['F']:.16f}",
        f"G_approx(T) (eV)           : {result['G_approx']:.16f}",
        "",
        "[Generated files]",
        "Initial structure          : initial_patch.xyz",
        "Relaxed structure          : relaxed_patch.xyz",
        "Relax log                  : relax_patch.log",
        "Relax trajectory           : relax_patch.traj",
        "Vibrational summary        : vibrations_summary.txt",
        "Temperature table          : patch_thermo_temperature.dat",
        "",
        "[Important note]",
        "This workflow uses HarmonicThermo, not IdealGasThermo.",
        "For supported fragments, fix boundary atoms or restrict vibrated atoms to avoid rigid-body modes.",
    ]

    write_output(output_dir / "patch_thermo_summary.out", output_lines)

    if args.clean_vib:
        vib.clean()

    print(f"Patch thermodynamic summary written to: {output_dir / 'patch_thermo_summary.out'}")
