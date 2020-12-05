"""Command line interface."""

import json
import click
import gpb.templating as templating


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


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
