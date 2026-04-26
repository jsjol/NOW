"""Save and load NOW_config to/from JSON and YAML."""
from __future__ import annotations

import json
from pathlib import Path

from ..config import NOW_config


def save_config(config: NOW_config, path: str | Path) -> None:
    """Save config to JSON or YAML (determined by file extension)."""
    path = Path(path)
    d = config.to_dict()

    if path.suffix in ('.yaml', '.yml'):
        try:
            import yaml
        except ImportError as e:
            raise ImportError("pyyaml is required for YAML format: pip install pyyaml") from e
        with open(path, 'w') as f:
            yaml.safe_dump(d, f, default_flow_style=False, sort_keys=False)
    else:
        with open(path, 'w') as f:
            json.dump(d, f, indent=2)


def load_config(path: str | Path) -> NOW_config:
    """Load config from JSON or YAML file."""
    path = Path(path)

    if path.suffix in ('.yaml', '.yml'):
        try:
            import yaml
        except ImportError as e:
            raise ImportError("pyyaml is required for YAML format: pip install pyyaml") from e
        with open(path) as f:
            d = yaml.safe_load(f)
    else:
        with open(path) as f:
            d = json.load(f)

    return NOW_config.from_dict(d)
