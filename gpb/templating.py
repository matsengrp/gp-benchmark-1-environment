"""Template input files and scripts."""

import math
import os
from jinja2 import Environment, FileSystemLoader, StrictUndefined


def template_file(template_name, template_dir, settings_dict, dest_path):
    """Build from a template."""
    env = Environment(loader=FileSystemLoader(template_dir), undefined=StrictUndefined)
    template = env.get_template(template_name)
    with open(dest_path, "w") as file_obj:
        file_obj.write(template.render(**settings_dict) + "\n")


def expand_mb_settings(settings_dict):
    """Expand settings in a settings dictionary to the full suite of things we
    want for MB."""

    if "burnin" not in settings_dict:
        settings_dict["burnin"] = math.floor(
            settings_dict["ngen"] * settings_dict["burnin_frac"]
        )

    if "burnin_samples" not in settings_dict:
        settings_dict["burnin_samples"] = math.ceil(
            settings_dict["burnin"] / settings_dict["samplefreq"]
        )


def make_paths_absolute(settings_dict):
    """Make the value of any key that ends with `_path` and absolute path."""

    for key, value in settings_dict.items():
        if key.endswith("_path"):
            settings_dict[key] = os.path.abspath(value)
