"""Command line interface."""

import json
import click
import click_config_file
import gpb.compare
import gpb.ourlibsbn as ourlibsbn
import gpb.templating as templating


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


# Entry point
@click.group(
    context_settings={"help_option_names": ["-h", "--help"]},
    invoke_without_command=True,
)
def cli():
    """Entrypoint."""
    pass  # pylint: disable=unnecessary-pass


@cli.command()
@click.argument("template_name", required=True)
@click.argument("settings_json", required=True, type=click.Path(exists=True))
@click.argument("dest_path", required=True, type=click.Path(exists=False))
def template(template_name, settings_json, dest_path):
    """Generate a file using one of our templates and the settings."""
    with open(settings_json, "r") as file:
        settings_dict = json.load(file)
    templating.template_file(template_name, settings_dict, dest_path)


@cli.command()
@click.argument("tree_path", required=True, type=click.Path(exists=True))
@click.argument("fasta_path", required=True, type=click.Path(exists=True))
@click.argument("out_csv_path", required=True, type=click.Path())
@click.option("--tol", type=float, default=1e-2)
@click.option("--max-iter", type=int, default=10)
@click_config_file.configuration_option(implicit=False, provider=json_provider)
def fit(tree_path, fasta_path, out_csv_path, tol, max_iter):
    """Fit an SBN using generalized pruning.

    Tree file is assumed to be a Newick file unless the path ends with `.nexus`.
    """
    ourlibsbn.gp_fit(tree_path, fasta_path, out_csv_path, tol, max_iter)


@cli.command()
@click.argument("newick_path", required=True, type=click.Path(exists=True))
@click.argument("out_csv_path", required=True, type=click.Path())
def sa(newick_path, out_csv_path):  # pylint: disable=invalid-name
    """Fit an SBN using simple average (SA) training.
    """
    ourlibsbn.simple_average(newick_path, out_csv_path)


@cli.command()
@click.argument("gp_csv", required=True, type=click.Path(exists=True))
@click.argument("sa_csv", required=True, type=click.Path(exists=True))
@click.argument("out_prefix", required=True, type=click.Path())
def compare(gp_csv, sa_csv, out_prefix):
    """Compare parameters between GP and SA.
    """
    gpb.compare.compare_parameters(gp_csv, sa_csv, out_prefix)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
