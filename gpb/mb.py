"""Template MrBayes input files."""

import math
import os
from jinja2 import Environment, PackageLoader, StrictUndefined


def generate_mb_input(settings_dict, dest_path):
    """
    Make MrBayes input file from a template.
    """
    env = Environment(
        loader=PackageLoader("gpb", "templates"), undefined=StrictUndefined
    )
    expand_mb_settings(settings_dict)
    make_paths_absolute(settings_dict)
    template = env.get_template(settings_dict["template"])
    with open(dest_path, "w") as file_obj:
        file_obj.write(template.render(**settings_dict) + "\n")


def expand_mb_settings(settings_dict):
    """
    Expand settings in a settings dictionary to the full suite of things expected by MB.
    """

    if "burnin" not in settings_dict:
        settings_dict["burnin"] = math.floor(
            settings_dict["ngen"] * settings_dict["burnin_frac"]
        )


def make_paths_absolute(settings_dict):
    """
    Make the value of any key that ends with `_path` and abosolute path.
    """

    for key, value in settings_dict.items():
        if key.endswith("_path"):
            settings_dict[key] = os.path.abspath(value)
