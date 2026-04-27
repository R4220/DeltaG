"""Helpers for configuration-file driven CLI execution.

This module lets the package keep a compact command line while moving most
workflow parameters into JSON or YAML configuration files.
"""

import argparse
import json
from pathlib import Path


CONFIG_COMMAND_KEYS = ("command", "workflow")
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

    return resolve_config_paths(data, path.parent)


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
