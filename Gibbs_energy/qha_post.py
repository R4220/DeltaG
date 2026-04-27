"""Post-processing tools for ``phonopy-qha`` output.

This module does not run phonons itself. Instead, it reconstructs convenient
thermodynamic tables starting from an existing ``phonopy-qha.out`` file and,
optionally, from a lightweight summary file with metadata such as pressure and
formula units.
"""

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .constants import PRESSURE_GPA_TO_EV_A3
from .io_utils import write_table


def parse_summary(path: Optional[Path]) -> Dict[str, str]:
    """Read a simple key-value summary file."""
    if path is None or not path.exists():
        return {}

    data: Dict[str, str] = {}

    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in line or line.lstrip().startswith("#"):
            continue

        key, value = line.split(":", 1)
        data[key.strip()] = value.strip()

    return data


def parse_phonopy_qha(path: Path) -> Tuple[List[float], ...]:
    """Read phonopy-qha.out columns: T, G, B0, B0prime, V."""
    rows = []

    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()

        if len(parts) != 5:
            continue

        try:
            rows.append([float(value) for value in parts])
        except ValueError:
            continue

    if not rows:
        raise ValueError(f"No QHA table rows found in {path}")

    temperatures = [row[0] for row in rows]
    gibbs = [row[1] for row in rows]
    bulk_modulus = [row[2] for row in rows]
    bulk_derivative = [row[3] for row in rows]
    volumes = [row[4] for row in rows]

    return temperatures, gibbs, bulk_modulus, bulk_derivative, volumes


def first_derivative(x_values: Sequence[float], y_values: Sequence[float]) -> List[float]:
    """Compute dy/dx on a discrete grid using finite differences."""
    if len(x_values) != len(y_values):
        raise ValueError("x_values and y_values must have the same length.")

    if len(x_values) < 2:
        raise ValueError("At least two points are required to compute a derivative.")

    derivatives: List[float] = []
    n_points = len(x_values)

    for index in range(n_points):
        # Use one-sided finite differences at the boundaries and a centered
        # formula in the interior.
        if index == 0:
            if n_points >= 3:
                dx = x_values[1] - x_values[0]
                value = (-3.0 * y_values[0] + 4.0 * y_values[1] - y_values[2]) / (2.0 * dx)
            else:
                value = (y_values[1] - y_values[0]) / (x_values[1] - x_values[0])

        elif index == n_points - 1:
            if n_points >= 3:
                dx = x_values[-1] - x_values[-2]
                value = (3.0 * y_values[-1] - 4.0 * y_values[-2] + y_values[-3]) / (2.0 * dx)
            else:
                value = (y_values[-1] - y_values[-2]) / (x_values[-1] - x_values[-2])

        else:
            value = (y_values[index + 1] - y_values[index - 1]) / (
                x_values[index + 1] - x_values[index - 1]
            )

        derivatives.append(value)

    return derivatives


def write_detailed_qha_table(
    path: Path,
    temperatures: Sequence[float],
    enthalpy: Sequence[float],
    gibbs: Sequence[float],
    entropy: Sequence[float],
    volumes: Sequence[float],
    bulk_modulus: Sequence[float],
    bulk_derivative: Sequence[float],
    pressure_gpa: Optional[float],
) -> None:
    """Write a detailed QHA thermodynamic table."""
    pressure_eV_A3 = None if pressure_gpa is None else pressure_gpa * PRESSURE_GPA_TO_EV_A3

    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# temperature[K] H[eV/cell] G[eV/cell] S[eV/K/cell] "
            "T*S[eV/cell] Veq[A^3/cell] B0[GPa] B0prime"
        )

        if pressure_eV_A3 is not None:
            handle.write(" pV[eV/cell] F[eV/cell]")

        handle.write("\n")
        handle.write("# S(T,p) is reconstructed as -dG/dT from the discrete QHA Gibbs curve.\n")

        for idx, temperature in enumerate(temperatures):
            ts_term = temperature * entropy[idx]

            line = (
                f"{temperature:16.6f} "
                f"{enthalpy[idx]:18.10f} "
                f"{gibbs[idx]:18.10f} "
                f"{entropy[idx]:18.10f} "
                f"{ts_term:18.10f} "
                f"{volumes[idx]:16.8f} "
                f"{bulk_modulus[idx]:12.6f} "
                f"{bulk_derivative[idx]:12.6f}"
            )

            if pressure_eV_A3 is not None:
                pv_term = pressure_eV_A3 * volumes[idx]
                helmholtz = gibbs[idx] - pv_term
                line += f" {pv_term:18.10f} {helmholtz:18.10f}"

            handle.write(line + "\n")


def run_qha_post(args) -> None:
    """Run QHA post-processing from phonopy-qha output."""
    phonopy_qha_path = Path(args.phonopy_qha).expanduser().resolve()
    summary_path = None if args.qha_summary is None else Path(args.qha_summary).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()

    output_dir.mkdir(parents=True, exist_ok=True)

    summary = parse_summary(summary_path)

    temperatures, gibbs, bulk_modulus, bulk_derivative, volumes = parse_phonopy_qha(
        phonopy_qha_path
    )

    # Thermodynamic entropy is reconstructed from the discrete Gibbs curve.
    dgdT = first_derivative(temperatures, gibbs)
    entropy = [-value for value in dgdT]

    enthalpy = [
        g_value + temperature * s_value
        for temperature, g_value, s_value in zip(temperatures, gibbs, entropy)
    ]

    pressure_gpa = None
    if "pressure_GPa" in summary:
        pressure_gpa = float(summary["pressure_GPa"])

    formula_units = None
    if "formula_units" in summary:
        formula_units = float(summary["formula_units"])

    write_table(
        output_dir / "enthalpy-temperature.dat",
        "# temperature[K] enthalpy[eV/cell]\n",
        temperatures,
        enthalpy,
    )

    write_table(
        output_dir / "entropy-temperature.dat",
        "# temperature[K] entropy[eV/K/cell]\n",
        temperatures,
        entropy,
    )

    write_table(
        output_dir / "gibbs-temperature.dat",
        "# temperature[K] gibbs[eV/cell]\n",
        temperatures,
        gibbs,
    )

    write_detailed_qha_table(
        output_dir / "qha_thermo_temperature.dat",
        temperatures,
        enthalpy,
        gibbs,
        entropy,
        volumes,
        bulk_modulus,
        bulk_derivative,
        pressure_gpa,
    )

    if formula_units is not None and formula_units != 0.0:
        write_table(
            output_dir / "enthalpy-temperature_per_formula.dat",
            "# temperature[K] enthalpy[eV/f.u.]\n",
            temperatures,
            [value / formula_units for value in enthalpy],
        )

        write_table(
            output_dir / "entropy-temperature_per_formula.dat",
            "# temperature[K] entropy[eV/K/f.u.]\n",
            temperatures,
            [value / formula_units for value in entropy],
        )

        write_table(
            output_dir / "gibbs-temperature_per_formula.dat",
            "# temperature[K] gibbs[eV/f.u.]\n",
            temperatures,
            [value / formula_units for value in gibbs],
        )

    print(f"Wrote QHA thermodynamic tables to: {output_dir}")
