"""Command line interface."""

import json
import click
import click_config_file
import gpb.compare
import gpb.outside
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
    """Top-level CLI for subcommands."""
    pass  # pylint: disable=unnecessary-pass


@cli.command()
@click.argument("template_name", required=True)
@click.argument("settings_json", required=True, type=click.Path(exists=True))
@click.argument("dest_path", required=True, type=click.Path(exists=False))
@click.option("--make-paths-absolute", is_flag=True, help="Make paths absolute.")
@click.option("--mb", is_flag=True, help="Perform MB-specific parameter expansion.")
def template(template_name, settings_json, dest_path, make_paths_absolute, mb):
    """Generate a file using one of our templates and the settings."""
    with open(settings_json, "r") as file:
        settings_dict = json.load(file)

    if make_paths_absolute:
        templating.expand_mb_settings(settings_dict)

    if mb:
        templating.make_paths_absolute(settings_dict)

    templating.template_file(template_name, settings_dict, dest_path)


@cli.command()
@click.argument("newick_path", required=True, type=click.Path(exists=True))
@click.argument("fasta_path", required=True, type=click.Path(exists=True))
@click.argument("out_csv_prefix", required=True, type=click.Path())
@click.option("--tol", type=float, default=1e-2)
@click.option("--max-iter", type=int, default=10)
@click.option("--mmap-path", type=click.Path(), default="mmap.dat")
@click_config_file.configuration_option(implicit=False, provider=json_provider)
def fit(newick_path, fasta_path, out_csv_prefix, tol, max_iter, mmap_path):
    """Fit an SBN using generalized pruning."""
    ourlibsbn.gp_fit(newick_path, fasta_path, out_csv_prefix, tol, max_iter, mmap_path)


@cli.command()
@click.argument("newick_path", required=True, type=click.Path(exists=True))
@click.argument("out_csv_prefix", required=True, type=click.Path())
def sa(newick_path, out_csv_prefix):  # pylint: disable=invalid-name
    """Fit an SBN using simple average (SA) training and spit the SBN parameters
    to CSV, as well as the unconditional subsplit probabilities to `.subsplit.csv`."""
    ourlibsbn.simple_average(newick_path, out_csv_prefix)


@cli.command()
@click.argument("gp_csv", required=True, type=click.Path(exists=True))
@click.argument("sa_csv", required=True, type=click.Path(exists=True))
@click.argument("sa_subsplit_csv", required=True, type=click.Path(exists=True))
@click.argument("out_prefix", required=True, type=click.Path())
def compare(gp_csv, sa_csv, sa_subsplit_csv, out_prefix):
    """Compare parameters between GP and SA, incorporating SA subsplit probabilities."""
    gpb.compare.compare_parameters(gp_csv, sa_csv, sa_subsplit_csv, out_prefix)


@cli.command()
@click.argument("newick_path", required=True, type=click.Path(exists=True))
@click.argument("sbn_parameter_csv", required=True, type=click.Path(exists=True))
@click.argument("out_csv_path", required=True, type=click.Path())
def treeprob(newick_path, sbn_parameter_csv, out_csv_path):
    """Calculate probabilities of the currently loaded trees and spit to CSV."""
    ourlibsbn.tree_probability(newick_path, sbn_parameter_csv, out_csv_path)


@cli.command()
@click.argument("newick_glob", required=True)
@click.argument("fasta_path", required=True, type=click.Path(exists=True))
@click.argument("out_csv_path", required=True, type=click.Path())
def treemarginal(newick_glob, fasta_path, out_csv_path):
    """Directly estimate the marginal log likelihood for trees supplied in a file for
    each file in the supplied Newick path glob."""
    ourlibsbn.tree_marginal(newick_glob, fasta_path, out_csv_path)


@cli.command()
@click.argument("newick_path", required=True, type=click.Path(exists=True))
@click.argument("fasta_path", required=True, type=click.Path(exists=True))
@click.argument("pcsp_csv_path", required=True, type=click.Path(exists=True))
@click.option("--tol", type=float, default=1e-2)
@click.option("--max-iter", type=int, default=10)
@click_config_file.configuration_option(
    implicit=False, cmd_name="fit", provider=json_provider
)
def treeexport(newick_path, fasta_path, pcsp_csv_path, tol, max_iter):
    """Export trees with GP branch lengths for each supplied PCSP."""
    ourlibsbn.export_trees_with_subsplits(
        newick_path, fasta_path, pcsp_csv_path, tol, max_iter
    )


@cli.command()
@click.argument("direct_marginals_csv", required=True, type=click.Path(exists=True))
@click.argument("prior_csv", required=True, type=click.Path(exists=True))
@click.argument("sa_csv", required=True, type=click.Path(exists=True))
@click.argument("out_prefix", required=True, type=click.Path())
def comparedirect(direct_marginals_csv, prior_csv, sa_csv, out_prefix):
    """Compare to the "direct" marginal likelihoods for estimating SBN parameters."""
    gpb.compare.compare_to_direct(direct_marginals_csv, prior_csv, sa_csv, out_prefix)


@cli.command()
@click.argument("original_path", required=True, type=click.Path(exists=True))
@click.argument("outside_path", required=True, type=click.Path(exists=True))
@click.argument("out_path", required=True, type=click.Path())
def outsideprob(original_path, outside_path, out_path):
    """Export the probability for the outside GPCSP with the largest parent subsplit."""
    gpb.outside.export_line_for_biggest_outside_gpcsp(
        original_path, outside_path, out_path
    )


@cli.command()
@click.argument("sbn_csv_path", required=True, type=click.Path(exists=True))
@click.argument("out_path", required=True, type=click.Path())
def addmeta(sbn_csv_path, out_path):
    """Add metadata about a subsplit csv and write it back out."""
    gpb.compare.add_metadata_to_sbn_csv(sbn_csv_path, out_path)


@cli.command()
@click.argument("outside_csv_path", required=True, type=click.Path(exists=True))
@click.argument("inside_csv_path", required=True, type=click.Path(exists=True))
@click.argument("out_path", required=True, type=click.Path())
def outsideplot(outside_csv_path, inside_csv_path, out_path):
    """Plot a comparison of outside and inside probabilities. We assume that all CSVs
    have "metadata" in the input CSV."""
    gpb.outside.plot_outside_prob_comparison(
        outside_csv_path, inside_csv_path, out_path
    )


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
