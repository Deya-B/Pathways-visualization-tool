import os
import yaml
import logging
from dataclasses import dataclass

@dataclass
class GPMLConfig:
    pathway_title: str
    organism: str
    id_data_file: str
    relations_file: str
    output_filename: str
    delimiter: str = "\t"
    logging_level: str = "INFO"
    logging_format: str = "%(levelname)s: %(message)s"


def load_config(path: str) -> GPMLConfig:
    """... {name} substitution implemented ..."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Config file not found: {path}")
    
    with open(path, "r", encoding="utf-8") as f:
        raw_cfg = yaml.safe_load(f) or {}

    # Validate required keys
    missing = []
    for key in ("pathway_title", "input", "output"):
        if key not in raw_cfg:
            missing.append(key)
    if missing:
        raise ValueError(f"Missing required config keys: {', '.join(missing)}")

    # variables available in templates
    vars_dict = {"name": raw_cfg.get("name", ""),}

    cfg = _substitute_vars(raw_cfg, vars_dict)

    in_cfg = cfg.get("input", {})
    out_cfg = cfg.get("output", {})
    log_cfg = cfg.get("logging", {})

    in_folder = in_cfg.get("folder", "")
    out_folder = out_cfg.get("folder", "")

    id_path = os.path.join(in_folder, in_cfg["id_data_file"])
    rel_path = os.path.join(in_folder, in_cfg["relations_file"])
    out_path = os.path.join(out_folder, out_cfg["filename"])

    return GPMLConfig(
        pathway_title=cfg["pathway_title"],
        organism=cfg.get("organism", "Homo sapiens"),
        id_data_file=id_path,
        relations_file=rel_path,
        output_filename=out_path,
        delimiter=in_cfg.get("delimiter", "\t"),
        logging_level=log_cfg.get("level", "INFO"),
        logging_format=log_cfg.get("format", "%(levelname)s: %(message)s"),
    )


def _substitute_vars(obj, vars_dict):
    if isinstance(obj, str):
        return obj.format(**vars_dict)
    if isinstance(obj, dict):
        return {k: _substitute_vars(v, vars_dict) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_substitute_vars(v, vars_dict) for v in obj]
    return obj


def setup_logging(config: GPMLConfig):
    level_name = (config.logging_level or "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    logging.basicConfig(level=level, format=config.logging_format)

