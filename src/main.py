"""Command-line entry point for the GPML pathway builder.

Parses command-line arguments, loads the YAML configuration, and runs
the GPML generation pipeline. Configuration and file-related errors
are logged and then re-raised.
"""

from pathways_tool.config import load_config
from pathways_tool.cli import run_from_config
import argparse
import logging

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build GPML pathway from ID and relations tables."
    )
    parser.add_argument("-c", "--config", required=True,
                        help="Path to YAML configuration file.")
    args = parser.parse_args()

    try:
        cfg = load_config(args.config)
        run_from_config(cfg)
    except FileNotFoundError as e:
        logging.error(str(e))
        raise
    except ValueError as e:
        logging.error(f"Configuration error: {e}")
        raise