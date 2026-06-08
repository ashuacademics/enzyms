import json
import os

import yaml


def load_config(filename):
    _, extension = os.path.splitext(filename.lower())
    with open(filename, "r") as file:
        if extension == ".json":
            return json.load(file)
        if extension in {".yaml", ".yml", ".param", ""}:
            return yaml.safe_load(file)
    raise ValueError(f"Unsupported parameter file extension: {extension}")
