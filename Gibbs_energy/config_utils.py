"""Helpers for configuration-file driven CLI execution.

This module lets the package keep a compact command line while moving most
workflow parameters into JSON or YAML configuration files.
"""

import argparse
import json
from pathlib import Path


CONFIG_COMMAND_KEYS = ("command",)
PATH_LIKE_KEYS = {
    "geometry_file",
    "model_path",
    "output_dir",
    "phonopy_qha",
    "qha_summary",
}


def load_config_file(config_path):
    """Load a JSON or YAML configuration file into a dictionary."""
    path = Path(config_path).expanduser().resolve()

    if not path.exists():
        raise ValueError(f"Configuration file not found: {path}")

    suffix = path.suffix.lower()
    text = path.read_text(encoding="utf-8")

    if suffix == ".json":
        data = json.loads(text)
    elif suffix in {".yaml", ".yml"}:
        try:
            import yaml
        except ImportError as exc:
            raise ValueError(
                "YAML support requires PyYAML. Install it or use a JSON config file."
            ) from exc

        data = yaml.safe_load(text)
    else:
        raise ValueError(
            f"Unsupported config format for {path}. Use .json, .yaml, or .yml."
        )

    if not isinstance(data, dict):
        raise ValueError(f"Configuration file must contain a top-level mapping: {path}")

    data = normalize_config_schema(data)
    return resolve_config_paths(data, path.parent)


def normalize_config_schema(config):
    """Normalize nested configs to the flat CLI-oriented structure."""
    config = dict(config)
    config.pop("schema_version", None)
    mode = config.get("mode")

    if mode is None:
        raise ValueError("Configuration files must define 'mode'.")

    if mode == "periodic":
        return normalize_periodic_config(config)

    if mode == "molecule":
        return normalize_molecule_config(config)

    if mode == "mixed":
        raise ValueError("Mode 'mixed' is planned but not implemented yet.")

    if mode in {"qha-post", "qha_post"}:
        return normalize_qha_post_config(config)

    raise ValueError(f"Unsupported mode '{mode}' in config file.")


def normalize_periodic_config(config):
    """Normalize the nested periodic config into flat CLI keys."""
    validate_allowed_keys(
        config,
        {"mode", "model", "thermo", "plots", "output", "periodic"},
        "top-level periodic config",
    )

    model = require_mapping(config, "model")
    periodic = require_mapping(config, "periodic")
    thermo = optional_mapping(config, "thermo")
    plots = optional_mapping(config, "plots")
    output = optional_mapping(config, "output")
    structure = require_mapping(periodic, "structure")
    relax = optional_mapping(periodic, "relax")
    phonons = optional_mapping(periodic, "phonons")
    temperature_grid = optional_mapping(periodic, "temperature_grid")

    validate_allowed_keys(model, {"path", "device"}, "model")
    validate_allowed_keys(thermo, {"temperature", "pressure"}, "thermo")
    validate_allowed_keys(
        plots,
        {"enabled", "phonon_dos", "band_structure", "thermo_curves", "format", "dpi"},
        "plots",
    )
    validate_allowed_keys(output, {"dir", "write_summary"}, "output")
    validate_allowed_keys(
        periodic,
        {"kind", "structure", "formula_units", "relax", "phonons", "temperature_grid"},
        "periodic",
    )
    validate_allowed_keys(structure, {"geometry_file", "lattice_geometry", "lattice_basis"}, "periodic.structure")
    validate_allowed_keys(relax, {"fmax"}, "periodic.relax")
    validate_allowed_keys(
        phonons,
        {"supercell", "delta", "dos_kpts", "bandpath", "emax"},
        "periodic.phonons",
    )
    validate_allowed_keys(
        temperature_grid,
        {"t_min", "t_max", "t_step"},
        "periodic.temperature_grid",
    )

    kind = periodic.get("kind")
    if kind not in {"bulk", "surface"}:
        raise ValueError("periodic.kind must be either 'bulk' or 'surface'.")

    if "geometry_file" in structure and (
        "lattice_geometry" in structure or "lattice_basis" in structure
    ):
        raise ValueError(
            "periodic.structure must define either geometry_file or lattice_geometry/lattice_basis, not both."
        )

    flat = {"command": "bulk"}

    copy_mapping(flat, model, {"path": "model_path", "device": "device"})
    copy_mapping(flat, thermo, {"temperature": "temperature", "pressure": "pressure"})
    copy_mapping(flat, output, {"dir": "output_dir"})
    copy_mapping(
        flat,
        periodic,
        {"formula_units": "formula_units"},
    )
    copy_mapping(flat, structure, {"geometry_file": "geometry_file"})

    if "lattice_geometry" in structure:
        flat["lattice_geometry"] = structure["lattice_geometry"]
    if "lattice_basis" in structure:
        flat["lattice_basis"] = structure["lattice_basis"]

    if ("lattice_geometry" in structure) != ("lattice_basis" in structure):
        raise ValueError(
            "periodic.structure must define lattice_geometry and lattice_basis together."
        )

    copy_mapping(flat, relax, {"fmax": "fmax"})
    copy_mapping(
        flat,
        phonons,
        {
            "supercell": "supercell",
            "delta": "delta",
            "dos_kpts": "dos_kpts",
            "bandpath": "bandpath",
            "emax": "emax",
        },
    )
    copy_mapping(
        flat,
        temperature_grid,
        {"t_min": "t_min", "t_max": "t_max", "t_step": "t_step"},
    )

    return flat


def normalize_molecule_config(config):
    """Normalize the nested molecular config into flat CLI keys."""
    validate_allowed_keys(
        config,
        {"mode", "model", "thermo", "plots", "output", "molecule"},
        "top-level molecule config",
    )

    model = require_mapping(config, "model")
    molecule = require_mapping(config, "molecule")
    thermo = optional_mapping(config, "thermo")
    plots = optional_mapping(config, "plots")
    output = optional_mapping(config, "output")
    relax = optional_mapping(molecule, "relax")
    vibrations = optional_mapping(molecule, "vibrations")

    validate_allowed_keys(model, {"path", "device"}, "model")
    validate_allowed_keys(thermo, {"temperature", "pressure"}, "thermo")
    validate_allowed_keys(
        plots,
        {"enabled", "phonon_dos", "band_structure", "thermo_curves", "format", "dpi"},
        "plots",
    )
    validate_allowed_keys(output, {"dir", "write_summary"}, "output")
    validate_allowed_keys(
        molecule,
        {
            "geometry_file",
            "mol_geometry",
            "symmetry_number",
            "spin",
            "vacuum",
            "relax",
            "vibrations",
        },
        "molecule",
    )
    validate_allowed_keys(relax, {"fmax"}, "molecule.relax")
    validate_allowed_keys(vibrations, {"clean"}, "molecule.vibrations")

    flat = {"command": "molecule"}

    copy_mapping(flat, model, {"path": "model_path", "device": "device"})
    copy_mapping(flat, thermo, {"temperature": "temperature", "pressure": "pressure"})
    copy_mapping(flat, output, {"dir": "output_dir"})
    copy_mapping(
        flat,
        molecule,
        {
            "geometry_file": "geometry_file",
            "mol_geometry": "mol_geometry",
            "symmetry_number": "symmetry_number",
            "spin": "spin",
            "vacuum": "vacuum",
        },
    )
    copy_mapping(flat, relax, {"fmax": "fmax"})

    if "clean" in vibrations:
        flat["clean_vib"] = vibrations["clean"]

    return flat


def normalize_qha_post_config(config):
    """Normalize the nested QHA post config into flat CLI keys."""
    validate_allowed_keys(
        config,
        {"mode", "output", "qha_post"},
        "top-level qha-post config",
    )

    output = optional_mapping(config, "output")
    qha_post = require_mapping(config, "qha_post")

    validate_allowed_keys(output, {"dir", "write_summary"}, "output")
    validate_allowed_keys(
        qha_post,
        {"phonopy_qha", "qha_summary"},
        "qha_post",
    )

    flat = {"command": "qha-post"}
    copy_mapping(flat, output, {"dir": "output_dir"})
    copy_mapping(
        flat,
        qha_post,
        {
            "phonopy_qha": "phonopy_qha",
            "qha_summary": "qha_summary",
        },
    )
    return flat


def require_mapping(mapping, key):
    """Return one required nested mapping from a config dictionary."""
    value = mapping.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Config block '{key}' must be a mapping.")
    return value


def optional_mapping(mapping, key):
    """Return one optional nested mapping from a config dictionary."""
    value = mapping.get(key)
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError(f"Config block '{key}' must be a mapping.")
    return value


def validate_allowed_keys(mapping, allowed_keys, label):
    """Validate that a mapping contains only known keys."""
    unknown_keys = sorted(set(mapping) - set(allowed_keys))
    if unknown_keys:
        unknown_list = ", ".join(unknown_keys)
        raise ValueError(f"Unknown key(s) in {label}: {unknown_list}")


def copy_mapping(destination, source, key_map):
    """Copy selected keys from one mapping into another."""
    for source_key, destination_key in key_map.items():
        if source_key in source:
            destination[destination_key] = source[source_key]


def resolve_config_paths(config, base_dir):
    """Resolve known path-like entries relative to the config file location."""
    resolved = dict(config)

    for key in PATH_LIKE_KEYS:
        value = resolved.get(key)

        if value is None or not isinstance(value, str):
            continue

        candidate = Path(value).expanduser()
        if not candidate.is_absolute():
            resolved[key] = str((base_dir / candidate).resolve())

    return resolved


def parse_override(override_text):
    """Parse one KEY=VALUE override from the command line."""
    if "=" not in override_text:
        raise ValueError(
            f"Invalid override '{override_text}'. Use the form KEY=VALUE."
        )

    key, raw_value = override_text.split("=", 1)
    key = key.strip()
    raw_value = raw_value.strip()

    if not key:
        raise ValueError(
            f"Invalid override '{override_text}'. The key cannot be empty."
        )

    try:
        value = json.loads(raw_value)
    except json.JSONDecodeError:
        value = raw_value

    return key, value


def apply_overrides(config, overrides):
    """Apply repeated KEY=VALUE command-line overrides to a config dict."""
    merged = dict(config)

    for override in overrides:
        key, value = parse_override(override)
        merged[key] = value

    return merged


def build_argv_from_config(parser, config):
    """Convert a config dictionary into a validated argparse-style argv list."""
    config = dict(config)
    command = config.get("command") or config.get("workflow")

    if command is None:
        raise ValueError(
            "Configuration file must define 'command' or 'workflow'."
        )

    subparser = get_subparser_for_command(parser, command)
    actions_by_dest = {
        action.dest: action
        for action in subparser._actions
        if action.option_strings and action.dest != "help"
    }

    allowed_keys = set(actions_by_dest) | set(CONFIG_COMMAND_KEYS)
    unknown_keys = sorted(set(config) - allowed_keys)

    if unknown_keys:
        unknown_list = ", ".join(unknown_keys)
        raise ValueError(
            f"Unknown configuration key(s) for command '{command}': {unknown_list}"
        )

    argv = [command]

    for key, value in config.items():
        if key in CONFIG_COMMAND_KEYS or value is None:
            continue

        action = actions_by_dest[key]
        argv.extend(action_tokens_from_value(action, value))

    return argv


def get_subparser_for_command(parser, command):
    """Return the argparse subparser corresponding to one workflow name."""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            if command not in action.choices:
                valid = ", ".join(sorted(action.choices))
                raise ValueError(
                    f"Unknown command '{command}' in config file. Valid commands: {valid}"
                )
            return action.choices[command]

    raise ValueError("Internal error: no workflow subparsers were defined.")


def action_tokens_from_value(action, value):
    """Serialize one config value into the command-line tokens for one action."""
    option = preferred_option_string(action)

    if isinstance(action, argparse._StoreTrueAction):
        if not isinstance(value, bool):
            raise ValueError(
                f"Config key '{action.dest}' must be true or false."
            )
        return [option] if value else []

    if action.nargs not in (None, "?"):
        if not isinstance(value, (list, tuple)):
            raise ValueError(
                f"Config key '{action.dest}' must be a list or tuple."
            )
        return [option] + [stringify_config_value(item) for item in value]

    return [option, stringify_config_value(value)]


def preferred_option_string(action):
    """Return the most readable option string for a parser action."""
    for option in action.option_strings:
        if option.startswith("--"):
            return option
    return action.option_strings[0]


def stringify_config_value(value):
    """Convert a Python config value back to a CLI token."""
    if isinstance(value, (dict, list)):
        return json.dumps(value)

    if isinstance(value, bool):
        return "true" if value else "false"

    return str(value)
