"""Resolve the configured solver output directory for an example run."""
from pathlib import Path
import re


def run_output_directory(run_directory):
    run_directory = Path(run_directory).resolve()
    input_path = run_directory / "IMEX_SfloW2D.inp"
    if not input_path.is_file():
        raise FileNotFoundError(f"solver input not found: {input_path}")

    parameters = {}
    for line in input_path.read_text(encoding="utf-8").splitlines():
        line = line.split("!", 1)[0]
        match = re.match(r"\s*(RUN_NAME|OUTPUT_BASE_DIRECTORY)\s*=\s*(.*?)\s*,?\s*$",
                         line, re.IGNORECASE)
        if match:
            parameters[match.group(1).upper()] = match.group(2).strip().strip("\"'")

    run_name = parameters.get("RUN_NAME")
    if not run_name:
        raise ValueError(f"RUN_NAME is missing from {input_path}")
    base_directory = Path(parameters.get("OUTPUT_BASE_DIRECTORY", "OUTPUT"))
    if not base_directory.is_absolute():
        base_directory = run_directory / base_directory
    return base_directory / run_name
