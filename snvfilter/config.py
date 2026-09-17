"""Configuration loading and merging module.

Priority: CLI args > user config > default config
"""

import os
import copy

import yaml


# Path to the package's default config
_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CONFIG_PATH = os.path.join(_PACKAGE_DIR, "..", "config", "default.yaml")
_DEFAULT_RULES_PATH = os.path.join(_PACKAGE_DIR, "..", "config", "rules.yaml")


def _deep_merge(base, override):
    """Recursively merge *override* into *base* and return *base*.

    Only dicts are merged recursively; everything else is overwritten.
    """
    for key, val in override.items():
        if key in base and isinstance(base[key], dict) and isinstance(val, dict):
            _deep_merge(base[key], val)
        else:
            base[key] = copy.deepcopy(val)
    return base


def load_yaml(path):
    """Load a YAML file and return the parsed dict."""
    if not path or not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    return data or {}


def load_config(config_file=None, cli_args=None):
    """Load and merge configuration layers.

    Parameters
    ----------
    config_file : str or None
        Path to a user YAML config file.
    cli_args : dict or None
        Dict of CLI argument overrides (flat or nested).

    Returns
    -------
    dict
        Merged configuration dict.
    """
    # Layer 1: defaults
    config = load_yaml(_DEFAULT_CONFIG_PATH)

    # Layer 2: user config
    if config_file:
        user_cfg = load_yaml(config_file)
        _deep_merge(config, user_cfg)

    # Layer 3: CLI overrides
    if cli_args:
        _deep_merge(config, cli_args)

    return config


def load_rules(rules_file=None):
    """Load filtering rules.

    Parameters
    ----------
    rules_file : str or None
        Path to a custom rules YAML. Falls back to the bundled default.

    Returns
    -------
    dict
    """
    path = rules_file if rules_file and os.path.exists(rules_file) else _DEFAULT_RULES_PATH
    return load_yaml(path)
