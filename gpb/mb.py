"""Template MrBayes input files."""

from jinja2 import Environment, PackageLoader, StrictUndefined


def generate_mb_input(settings_dict, dest_path):
    """
    Make MrBayes input file from a template.
    """
    env = Environment(
        loader=PackageLoader('gpb', 'templates'),
        undefined=StrictUndefined
    )
    template = env.get_template('basic.mb')
    with open(dest_path, "w") as file_obj:
        file_obj.write(template.render(**settings_dict) + "\n")
