from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import tomllib
import yaml


class ConfigError(ValueError):
    pass


def load_config(path: str | Path) -> dict[str, Any]:
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise ConfigError(f"Config file not found: {p}")
    suffix = p.suffix.lower()
    text = p.read_text()
    if suffix in {".yaml", ".yml"}:
        data = yaml.safe_load(text)
    elif suffix == ".toml":
        data = tomllib.loads(text)
    else:
        raise ConfigError(f"Unsupported config format for {p}. Use YAML or TOML.")
    if not isinstance(data, dict):
        raise ConfigError(f"Config root must be a mapping in {p}")
    data["__config_path__"] = str(p)
    return data


def add_common_overrides(ap: argparse.ArgumentParser) -> argparse.ArgumentParser:
    ap.add_argument("--config", default="configs/atlas_v5.yaml")
    ap.add_argument("--root", dest="root_override", default=None)
    ap.add_argument("--in", dest="in_override", default=None)
    ap.add_argument("--out", dest="out_override", default=None)
    ap.add_argument("--out-tag", dest="out_tag_override", default=None)
    return ap


def merge_cli_overrides(cfg: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    c = copy.deepcopy(cfg)
    c.setdefault("paths", {})
    if args.root_override:
        c["root"] = args.root_override
    if args.in_override:
        c["paths"]["input"] = args.in_override
    if args.out_override:
        c["paths"]["output"] = args.out_override
    if args.out_tag_override:
        c["out_tag"] = args.out_tag_override
    return c


def _root(cfg: dict[str, Any]) -> Path:
    base = cfg.get("root") or "."
    return Path(base).expanduser().resolve()


def resolve_path(cfg: dict[str, Any], value: str | Path) -> Path:
    p = Path(value).expanduser()
    return p if p.is_absolute() else (_root(cfg) / p).resolve()


def require_file(path: Path, label: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing required {label}: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"{label} exists but is not a file: {path}")
    return path


def stamp_provenance(sidecar: Path, cfg: dict[str, Any], resolved: dict[str, Any]) -> None:
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "config_path": cfg.get("__config_path__"),
        "root": str(_root(cfg)),
        "resolved": resolved,
    }
    sidecar.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
