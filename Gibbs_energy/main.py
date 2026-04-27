"""Command-line interface for the thermochemistry mini-toolkit.

The current CLI exposes four independent workflows:

- ``bulk`` for harmonic bulk thermodynamics at fixed volume/cell;
- ``patch`` for harmonic thermodynamics of finite fragments in vacuum;
- ``molecule`` for ideal-gas thermochemistry of isolated molecules;
- ``qha-post`` for post-processing external phonopy-QHA results.
"""

import argparse
import json

from .config_utils import apply_overrides, build_argv_from_config, load_config_file
from .patch_thermo import run_patch
from .qha_post import run_qha_post
from .thermo_bulk import run_bulk
from .molecule_thermo import run_molecule


def build_bootstrap_parser():
    """Build the small parser used to detect config-driven execution."""
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument(
        "--config",
        default=None,
        help="Path to a JSON/YAML config file using the nested workflow structure.",
    )
    parser.add_argument(
        "--override",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Override one config entry from the command line. Can be repeated.",
    )

    return parser


def add_bulk_parser(subparsers):
    """Add command-line options for the bulk workflow."""
    parser = subparsers.add_parser(
        "bulk",
        help="Run ASE/MACE bulk phonon thermochemistry.",
    )

    parser.add_argument("--model-path", required=True, help="Path to the MACE model.")
    parser.add_argument("--device", default="cuda", help="Device used by MACE: cuda or cpu.")

    parser.add_argument(
        "--geometry-file",
        default=None,
        help="Periodic structure file readable by ASE.",
    )

    parser.add_argument(
        "--lattice-geometry",
        type=json.loads,
        default=None,
        help="JSON with a,b,c,alpha,beta,gamma,Inum.",
    )

    parser.add_argument(
        "--lattice-basis",
        type=json.loads,
        default=None,
        help="JSON list of atomic symbols and fractional positions.",
    )

    parser.add_argument(
        "--formula-units",
        type=int,
        default=None,
        help="Formula units in the simulation cell. If omitted, inferred automatically.",
    )
    parser.add_argument(
        "--kind",
        default="bulk",
        choices=["bulk", "surface"],
        help="Periodic system kind: bulk or surface.",
    )

    parser.add_argument("--fmax", type=float, default=0.01, help="Relaxation force threshold eV/A.")
    parser.add_argument("--supercell", type=int, nargs=3, default=(5, 5, 5))
    parser.add_argument("--delta", type=float, default=0.05, help="Finite displacement in A.")
    parser.add_argument("--dos-kpts", type=int, nargs=3, default=(40, 40, 40))
    parser.add_argument("--bandpath", default=None, help="ASE bandpath string, e.g. GMKG.")
    parser.add_argument("--emax", type=float, default=0.035, help="Max phonon energy in plot.")

    parser.add_argument("--temperature", type=float, default=298.15)
    parser.add_argument("--pressure", type=float, default=0.0, help="Pressure in GPa.")

    parser.add_argument("--t-min", type=float, default=0.0)
    parser.add_argument("--t-max", type=float, default=1000.0)
    parser.add_argument("--t-step", type=float, default=50.0)

    parser.add_argument("--output-dir", default="bulk_results")
    parser.add_argument(
        "--qha",
        action="store_true",
        help="Also run QHA post-processing at the end using an existing phonopy-qha.out file.",
    )
    parser.add_argument(
        "--phonopy-qha",
        default=None,
        help="Path to an existing phonopy-qha.out file used when --qha is enabled.",
    )
    parser.add_argument(
        "--qha-summary",
        default=None,
        help="Optional qha_summary.out used when --qha is enabled. If omitted, one is generated automatically.",
    )
    parser.add_argument(
        "--qha-output-dir",
        default=None,
        help="Optional output directory for QHA tables. Defaults to <bulk output>/qha_tables.",
    )

    parser.set_defaults(func=run_bulk)


def add_qha_post_parser(subparsers):
    """Add command-line options for the QHA post-processing workflow."""
    parser = subparsers.add_parser(
        "qha-post",
        help="Post-process phonopy-qha.out into thermodynamic tables.",
    )

    parser.add_argument(
        "--phonopy-qha",
        required=True,
        help="Path to phonopy-qha.out.",
    )

    parser.add_argument(
        "--qha-summary",
        default=None,
        help="Optional qha_summary.out with pressure_GPa and formula_units.",
    )

    parser.add_argument(
        "--output-dir",
        default="qha_tables",
        help="Directory where QHA tables are written.",
    )

    parser.set_defaults(func=run_qha_post)


def add_patch_parser(subparsers):
    """Add command-line options for harmonic patch thermochemistry."""
    parser = subparsers.add_parser(
        "patch",
        help="Run harmonic thermochemistry for a finite patch in vacuum.",
    )

    parser.add_argument("--model-path", required=True, help="Path to the MACE model.")
    parser.add_argument("--device", default="cuda", help="Device used by MACE: cuda or cpu.")

    parser.add_argument(
        "--geometry-file",
        required=True,
        help="Finite structure file readable by ASE.",
    )

    parser.add_argument(
        "--vacuum",
        type=float,
        default=15.0,
        help="Vacuum padding added around the finite structure in Angstrom.",
    )

    parser.add_argument(
        "--fixed-indices",
        type=int,
        nargs="+",
        default=None,
        help="Optional 0-based atom indices kept fixed during relaxation and vibrations.",
    )

    parser.add_argument(
        "--vibration-indices",
        type=int,
        nargs="+",
        default=None,
        help="Optional 0-based atom indices included in the finite-displacement vibrational analysis.",
    )

    parser.add_argument("--fmax", type=float, default=0.01, help="Relaxation force threshold eV/A.")
    parser.add_argument("--delta", type=float, default=0.01, help="Finite displacement in A for vibrations.")

    parser.add_argument("--temperature", type=float, default=298.15)
    parser.add_argument("--t-min", type=float, default=0.0)
    parser.add_argument("--t-max", type=float, default=1000.0)
    parser.add_argument("--t-step", type=float, default=50.0)

    parser.add_argument("--output-dir", default="patch_results")

    parser.add_argument(
        "--clean-vib",
        action="store_true",
        help="Remove ASE vibration displacement files after calculation.",
    )

    parser.set_defaults(func=run_patch)


def add_molecule_parser(subparsers):
    """Add command-line options for ideal-gas molecular thermochemistry."""
    parser = subparsers.add_parser(
        "molecule",
        help="Run ideal-gas molecular thermochemistry.",
    )

    parser.add_argument("--model-path", required=True, help="Path to the MACE model.")
    parser.add_argument("--device", default="cuda", help="Device used by MACE: cuda or cpu.")

    parser.add_argument(
        "--geometry-file",
        required=True,
        help="Molecular structure file readable by ASE.",
    )

    parser.add_argument(
        "--mol-geometry",
        required=True,
        choices=["monatomic", "linear", "nonlinear"],
        help="Molecular geometry type for IdealGasThermo.",
    )

    parser.add_argument(
        "--symmetry-number",
        type=int,
        required=True,
        help="Rotational symmetry number of the molecule.",
    )

    parser.add_argument(
        "--spin",
        type=float,
        default=0.0,
        help="Total electronic spin. Example: 1.0 for triplet O2, 0.0 for singlets.",
    )

    parser.add_argument("--vacuum", type=float, default=15.0)
    parser.add_argument("--fmax", type=float, default=0.01)

    parser.add_argument("--temperature", type=float, default=298.15)
    parser.add_argument(
        "--pressure",
        type=float,
        default=101325.0,
        help="Gas pressure in Pa.",
    )

    parser.add_argument("--output-dir", default="molecule_results")

    parser.add_argument(
        "--clean-vib",
        action="store_true",
        help="Remove ASE vibration displacement files after calculation.",
    )

    parser.set_defaults(func=run_molecule)

def build_workflow_parser(bootstrap_parser):
    """Build the main workflow parser used by both direct CLI and config mode."""
    parser = argparse.ArgumentParser(
        parents=[bootstrap_parser],
        description="Thermochemistry workflows for bulk crystals, finite patches, molecules, and QHA post-processing."
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Each subparser defines one scientific workflow with its own parameters.
    add_bulk_parser(subparsers)
    add_patch_parser(subparsers)
    add_qha_post_parser(subparsers)
    add_molecule_parser(subparsers)

    return parser


def parse_args(argv=None):
    """Parse command-line arguments and return the selected workflow config."""
    bootstrap_parser = build_bootstrap_parser()
    bootstrap_args, remaining = bootstrap_parser.parse_known_args(argv)

    parser = build_workflow_parser(bootstrap_parser)

    if bootstrap_args.config is not None:
        if remaining:
            parser.error(
                "When --config is used, do not pass workflow arguments directly. "
                "Use --override KEY=VALUE for small command-line changes."
            )

        try:
            config_data = load_config_file(bootstrap_args.config)
            merged_config = apply_overrides(config_data, bootstrap_args.override)
            config_argv = build_argv_from_config(parser, merged_config)
        except ValueError as exc:
            parser.error(str(exc))

        args = parser.parse_args(config_argv)
        args.config_file = bootstrap_args.config
        args.config_overrides = list(bootstrap_args.override)
        return args

    if bootstrap_args.override:
        parser.error("--override can only be used together with --config.")

    return parser.parse_args(argv)


def main():
    """Run the workflow selected on the command line."""
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
