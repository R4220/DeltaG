"""Bulk harmonic thermodynamics for periodic systems.

This module handles the bulk workflow in the harmonic approximation at fixed
cell shape and volume:

1. build or load a periodic structure;
2. relax atomic positions;
3. compute phonons with finite displacements;
4. derive thermodynamic quantities from the phonon DOS.

It does not perform thermal expansion. For that, the QHA workflow should be
used as a later step.
"""

from argparse import Namespace
from pathlib import Path

import numpy as np
from ase import units
from ase.io import write
from ase.optimize import QuasiNewton
from ase.thermochemistry import CrystalThermo

from .calculators import setup_calculator
from .constants import PRESSURE_GPA_TO_EV_A3
from .io_utils import write_output
from .phonons_bulk import plot_phonon_bs_and_dos, run_phonons
from .qha_post import run_qha_post
from .structures import build_structure, formula_strings, infer_formula_units


def crystal_energy_breakdown(phonon_energies, phonon_dos, potentialenergy, formula_units, temperature):
    """Compute E_pot, ZPE, thermal vibrational energy and U per formula unit."""
    omega_e = np.asarray(phonon_energies, dtype=float)
    dos_e = np.asarray(phonon_dos, dtype=float)

    # The first DOS point can be exactly zero because of the acoustic mode.
    if len(omega_e) > 0 and omega_e[0] == 0.0:
        omega_e = omega_e[1:]
        dos_e = dos_e[1:]

    e_pot = potentialenergy / formula_units

    e_zpe = np.trapz(0.5 * omega_e * dos_e, omega_e) / formula_units

    if temperature <= 0.0:
        e_thermal_vib = 0.0
    else:
        beta = 1.0 / (units.kB * temperature)
        vib_density = omega_e / (np.exp(omega_e * beta) - 1.0)
        e_thermal_vib = np.trapz(vib_density * dos_e, omega_e) / formula_units

    return {
        "E_pot": e_pot,
        "E_ZPE": e_zpe,
        "E_thermal_vib": e_thermal_vib,
        "U": e_pot + e_zpe + e_thermal_vib,
    }


def compute_bulk_thermo_at_temperature(
    thermo,
    phonon_energies,
    phonon_dos,
    potentialenergy,
    formula_units,
    temperature,
    pressure_gpa,
    volume_total,
):
    """Compute U, S, F, H and G at one temperature for the bulk workflow."""
    energy_terms = crystal_energy_breakdown(
        phonon_energies=phonon_energies,
        phonon_dos=phonon_dos,
        potentialenergy=potentialenergy,
        formula_units=formula_units,
        temperature=temperature,
    )

    U = energy_terms["U"]
    S = thermo.get_entropy(temperature=temperature, verbose=False)

    # Convert the user pressure to ASE-compatible energy-density units.
    pressure_eV_A3 = pressure_gpa * PRESSURE_GPA_TO_EV_A3
    pv_total = pressure_eV_A3 * volume_total
    pv_per_formula = pv_total / formula_units

    TS = temperature * S
    F = U - TS
    H = U + pv_per_formula
    G = F + pv_per_formula

    return {
        **energy_terms,
        "S": S,
        "TS": TS,
        "pV": pv_per_formula,
        "F": F,
        "H": H,
        "G": G,
    }


def write_bulk_temperature_table(
    output_path,
    temperatures,
    thermo,
    phonon_energies,
    phonon_dos,
    potentialenergy,
    formula_units,
    pressure_gpa,
    volume_total,
):
    """Write bulk U, S, F, H and G as functions of temperature."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# T[K] E_pot[eV/f.u.] E_ZPE[eV/f.u.] "
            "E_thermal_vib[eV/f.u.] U[eV/f.u.] S[eV/K/f.u.] "
            "TS[eV/f.u.] pV[eV/f.u.] F[eV/f.u.] H[eV/f.u.] G[eV/f.u.]\n"
        )

        for temperature in temperatures:
            result = compute_bulk_thermo_at_temperature(
                thermo=thermo,
                phonon_energies=phonon_energies,
                phonon_dos=phonon_dos,
                potentialenergy=potentialenergy,
                formula_units=formula_units,
                temperature=temperature,
                pressure_gpa=pressure_gpa,
                volume_total=volume_total,
            )

            handle.write(
                f"{temperature:16.6f} "
                f"{result['E_pot']:18.10f} "
                f"{result['E_ZPE']:18.10f} "
                f"{result['E_thermal_vib']:18.10f} "
                f"{result['U']:18.10f} "
                f"{result['S']:18.10f} "
                f"{result['TS']:18.10f} "
                f"{result['pV']:18.10f} "
                f"{result['F']:18.10f} "
                f"{result['H']:18.10f} "
                f"{result['G']:18.10f}\n"
            )


def maybe_run_qha_post(args, output_dir, formula_units):
    """Optionally run QHA post-processing after the bulk workflow."""
    if not getattr(args, "qha", False):
        return None, None

    if args.phonopy_qha is None:
        raise ValueError(
            "Bulk workflow with qha enabled requires phonopy-qha.out via qha.phonopy_qha or --phonopy-qha."
        )

    qha_output_dir = output_dir / "qha_tables"
    if args.qha_output_dir is not None:
        qha_output_dir = Path(args.qha_output_dir)

    qha_summary_path = args.qha_summary
    if qha_summary_path is not None:
        qha_summary_path = Path(qha_summary_path)

    if qha_summary_path is None or not qha_summary_path.exists():
        if qha_summary_path is None:
            qha_summary_path = output_dir / "qha_summary.out"
        write_output(
            qha_summary_path,
            [
                f"pressure_GPa: {args.pressure}",
                f"formula_units: {formula_units}",
            ],
        )

    qha_args = Namespace(
        phonopy_qha=args.phonopy_qha,
        qha_summary=str(qha_summary_path),
        output_dir=str(qha_output_dir),
    )
    run_qha_post(qha_args)
    return qha_output_dir, Path(qha_args.qha_summary)


def run_bulk(args):
    """Run the complete bulk ASE/MACE phonon thermochemistry workflow."""
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    calc = setup_calculator(args.model_path, device=args.device)
    atoms = build_structure(args, calc)

    # Normalize all extensive quantities to a reduced formula unit when possible.
    formula_units = args.formula_units or infer_formula_units(atoms)

    write(output_dir / "initial_crystal.xyz", atoms)

    optimizer = QuasiNewton(
        atoms,
        logfile=str(output_dir / "relax.log"),
        trajectory=str(output_dir / "relax.traj"),
    )
    optimizer.run(fmax=args.fmax)

    write(output_dir / "relaxed_crystal.xyz", atoms)

    potentialenergy = atoms.get_potential_energy()
    volume_total = atoms.get_volume()

    # Phonons are evaluated on the relaxed geometry at fixed cell.
    _, dos, band_structure = run_phonons(
        atoms=atoms,
        calculator=calc,
        supercell=args.supercell,
        delta=args.delta,
        dos_kpts=args.dos_kpts,
        bandpath=args.bandpath,
        output_dir=output_dir,
    )

    phonon_energies = dos.get_energies()
    phonon_dos = dos.get_weights()

    plot_phonon_bs_and_dos(
        dos=dos,
        band_structure=band_structure,
        output_path=output_dir / "phonon_BS_and_DOS.png",
        emax=args.emax,
    )

    thermo = CrystalThermo(
        phonon_energies=phonon_energies,
        phonon_DOS=phonon_dos,
        potentialenergy=potentialenergy,
        formula_units=formula_units,
    )

    temperatures = np.arange(args.t_min, args.t_max + args.t_step, args.t_step)

    write_bulk_temperature_table(
        output_path=output_dir / "bulk_thermo_temperature.dat",
        temperatures=temperatures,
        thermo=thermo,
        phonon_energies=phonon_energies,
        phonon_dos=phonon_dos,
        potentialenergy=potentialenergy,
        formula_units=formula_units,
        pressure_gpa=args.pressure,
        volume_total=volume_total,
    )

    result = compute_bulk_thermo_at_temperature(
        thermo=thermo,
        phonon_energies=phonon_energies,
        phonon_dos=phonon_dos,
        potentialenergy=potentialenergy,
        formula_units=formula_units,
        temperature=args.temperature,
        pressure_gpa=args.pressure,
        volume_total=volume_total,
    )

    qha_output_dir, qha_summary_path = maybe_run_qha_post(
        args=args,
        output_dir=output_dir,
        formula_units=formula_units,
    )

    cell_formula, reduced_formula = formula_strings(atoms, formula_units)

    system_lines = [
        f"Geometry file              : {args.geometry_file}",
        f"Relaxed cell formula       : {cell_formula}",
        f"Reduced formula unit       : {reduced_formula}",
        f"Formula units in cell      : {formula_units}",
    ]

    output_lines = [
        "# Bulk harmonic thermochemistry report",
        "",
        "[System]",
        *system_lines,
        "",
        "[Settings]",
        f"Model path                 : {args.model_path}",
        f"Device                     : {args.device}",
        f"Temperature selected (K)   : {args.temperature}",
        f"Pressure (GPa)             : {args.pressure}",
        f"Force threshold (eV/A)     : {args.fmax}",
        f"Phonon supercell           : {tuple(args.supercell)}",
        f"Displacement delta (A)     : {args.delta}",
        f"DOS k-point mesh           : {tuple(args.dos_kpts)}",
        f"Band path                  : {args.bandpath}",
        "",
        "[Results per formula unit]",
        f"Volume (A^3/f.u.)          : {volume_total / formula_units:.16f}",
        f"E_pot (eV/f.u.)            : {result['E_pot']:.16f}",
        f"E_ZPE (eV/f.u.)            : {result['E_ZPE']:.16f}",
        f"E_thermal_vib (eV/f.u.)    : {result['E_thermal_vib']:.16f}",
        f"U(T) (eV/f.u.)             : {result['U']:.16f}",
        f"S(T) (eV/K/f.u.)           : {result['S']:.16f}",
        f"T*S (eV/f.u.)              : {result['TS']:.16f}",
        f"pV (eV/f.u.)               : {result['pV']:.16f}",
        f"F(T) (eV/f.u.)             : {result['F']:.16f}",
        f"H(T,p) (eV/f.u.)           : {result['H']:.16f}",
        f"G(T,p) (eV/f.u.)           : {result['G']:.16f}",
        "",
        "[Generated files]",
        "Initial structure          : initial_crystal.xyz",
        "Relaxed structure          : relaxed_crystal.xyz",
        "Relax log                  : relax.log",
        "Relax trajectory           : relax.traj",
        "Phonon plot                : phonon_BS_and_DOS.png",
        "Temperature table          : bulk_thermo_temperature.dat",
    ]

    if qha_output_dir is not None:
        output_lines.extend(
            [
                f"QHA summary metadata       : {qha_summary_path}",
                f"QHA tables directory       : {qha_output_dir}",
            ]
        )

    output_lines.extend(
        [
            "",
            "[Important note]",
            "This is a bulk harmonic calculation at fixed cell. For thermal expansion and true G(T,p), use QHA.",
        ]
    )

    write_output(output_dir / "bulk_summary.out", output_lines)

    print(f"Bulk thermodynamic summary written to: {output_dir / 'bulk_summary.out'}")
