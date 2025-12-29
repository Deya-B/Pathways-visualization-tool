import os
import yaml
import logging
from dataclasses import dataclass

@dataclass
class GPMLConfig:
    """Configuration for building a GPML pathway.

    Attributes
    ----------
    pathway_title : str
        Title to use for the GPML pathway.
    organism : str
        Organism name to record in the GPML metadata.
    id_data_file : str
        Path to the TSV file with ID metadata (metabolites, pathways, enzymes).
    relations_file : str
        Path to the TSV file describing source–target relationships.
    output_filename : str
        Path where the resulting GPML file will be written.
    delimiter : str, optional
        Field delimiter used in the input TSV files, by default "\\t".
    logging_level : str, optional
        Logging level name (e.g., "INFO", "DEBUG"), by default "INFO".
    logging_format : str, optional
        Logging format string passed to logging.basicConfig, by default
        "%(levelname)s: %(message)s".
    """
    pathway_title: str
    organism: str
    id_data_file: str
    relations_file: str
    output_filename: str
    delimiter: str = "\t"
    logging_level: str = "INFO"
    logging_format: str = "%(levelname)s: %(message)s"


def load_config(path: str) -> GPMLConfig:
    """Load YAML configuration, apply {name} templates, and build GPMLConfig.

    The raw YAML configuration may contain `{name}` placeholders in string
    values. These are expanded using the required top-level ``name`` field
    before input/output paths and logging options are resolved. The function
    validates that mandatory top-level and nested keys are present and that
    ``name`` is non-empty.

    Parameters
    ----------
    path : str
        Path to a YAML configuration file describing input, output, and
        logging settings.

    Returns
    -------
    GPMLConfig
        Normalized configuration with resolved input/output file paths,
        logging options, and any `{name}` templates substituted.

    Raises
    ------
    FileNotFoundError
        If the configuration file does not exist.
    ValueError
        If required top-level keys (\"pathway_title\", \"name\", \"input\",
        \"output\") are missing, if ``name`` is empty, or if required keys
        are missing inside the ``input`` (folder, id_data_file,
        relations_file) or ``output`` (folder, filename) sections.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Config file not found: {path}")
    
    with open(path, "r", encoding="utf-8") as f:
        raw_cfg = yaml.safe_load(f) or {}

    # Validate required keys
    missing = []
    for key in ("pathway_title", "name", "input", "output"):
        if key not in raw_cfg:
            missing.append(key)
    if missing:
        raise ValueError (
            "Invalid config: missing top-level keys: "
            f"{', '.join(missing)}. "
            "Expected keys: pathway_title, input, name, "
            "output (and optional organism, logging)."
        )
    
    if not str(raw_cfg.get("name", "")).strip():
        raise ValueError(
            "Invalid config: 'name' is empty. "
            "This is used to build file names (e.g. {name}.tsv). "
            "Please set a non-empty 'name' value in the YAML."
        )

    # variables available in templates
    vars_dict = {"name": raw_cfg.get("name", ""),}

    cfg = _substitute_vars(raw_cfg, vars_dict)

    in_cfg = cfg.get("input", {})
    out_cfg = cfg.get("output", {})
    log_cfg = cfg.get("logging", {})

    # Nested keys under input/output
    missing_input = [k for k in ("folder", "id_data_file", "relations_file") if k not in in_cfg]
    missing_output = [k for k in ("folder", "filename") if k not in out_cfg]

    errors = []
    if missing_input:
        errors.append(
            "input section is missing: " + ", ".join(missing_input)
            + " (required: folder, id_data_file, relations_file)"
        )
    if missing_output:
        errors.append("output section is missing: " + ", ".join(missing_output)
                      + " (required: folder, filename)")
    if errors:
        raise ValueError("Invalid config:\n  - " + "\n  - ".join(errors))

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
    """Recursively substitute template variables in a nested config object.

    Parameters
    ----------
    obj : Any
        Configuration structure (string, dict, list, or other) potentially
        containing `{name}`-style format fields.
    vars_dict : dict
        Mapping of placeholder names to replacement values used with
        `str.format`.

    Returns
    -------
    Any
        Object of the same shape as `obj` with string values formatted
        using `vars_dict`.
    """
    if isinstance(obj, str):
        return obj.format(**vars_dict)
    if isinstance(obj, dict):
        return {k: _substitute_vars(v, vars_dict) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_substitute_vars(v, vars_dict) for v in obj]
    return obj


def setup_logging(config: GPMLConfig):
    """Configure the root logger from a GPMLConfig instance.

    Parameters
    ----------
    config : GPMLConfig
        Configuration object providing the logging level name and
        logging format string used to initialize `logging.basicConfig`.
    """
    level_name = (config.logging_level or "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    logging.basicConfig(level=level, format=config.logging_format)