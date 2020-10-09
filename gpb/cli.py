"""Command line interface."""

import pathlib
import json
import click
import click_config_file
import gpb.mb as mb

def alignment_path_of_prefix(prefix):
    return f"{prefix}.fasta"


def json_provider(file_path, cmd_name):
    """Enable loading of flags from a JSON file via click_config_file."""
    if cmd_name:
        with open(file_path) as config_data:
            config_dict = json.load(config_data)
            if cmd_name not in config_dict:
                if "default" in config_dict:
                    return config_dict["default"]
                # else:
                raise IOError(
                    f"Could not find a '{cmd_name}' or 'default' section in '{file_path}'"
                )
            return config_dict[cmd_name]
    # else:
    return None


def print_method_name_and_locals(method_name, local_variables):
    """Print method name and local variables."""
    if "ctx" in local_variables:
        del local_variables["ctx"]
    print(f"{method_name}{local_variables})")


def restrict_dict_to_params(d_to_restrict, cmd):
    """Restrict the given dictionary to the names of parameters for cmd."""
    param_names = {param.name for param in cmd.params}
    return {key: d_to_restrict[key] for key in d_to_restrict if key in param_names}


def dry_run_option(command):
    return click.option(
        "--dry-run",
        is_flag=True,
        help="Only print paths and files to be made, rather than actually making them.",
    )(command)


# Entry point
@click.group(
    context_settings={"help_option_names": ["-h", "--help"]},
    invoke_without_command=True,
)
def cli():
    pass


@cli.command()
@click.argument("dest_path", required=True, type=click.Path(exists=False))
def template(dest_path):
    """Generate a MrBayes input file."""
    settings_dict = {"ngen":1234, "burnin":3}
    mb.generate_mb_input(settings_dict, dest_path)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
