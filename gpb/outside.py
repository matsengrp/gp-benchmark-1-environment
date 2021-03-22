"""Functions to compare subsplits outside the support."""

import click
import pandas as pd


def parent_split_size(gpcsp):
    "Count the number of taxa in the parent split of the PCSP."
    gpcsp_list = gpcsp.split("|")
    assert len(gpcsp_list) == 3
    return (gpcsp_list[0] + gpcsp_list[1]).count("1")


def get_biggest_gpcsp(gpcsp_list):
    "Return the GPCSP with the most taxa in its parent split."
    split_size_df = (
        pd.DataFrame(
            {
                "gpcsp": gpcsp_list,
                "parent_split_size": [parent_split_size(gpcsp) for gpcsp in gpcsp_list],
            }
        )
        .sort_values("parent_split_size", ascending=False)
        .reset_index()
    )
    split_sizes = split_size_df["parent_split_size"]
    assert len(split_sizes) > 0
    if len(split_sizes) > 1:
        assert split_sizes[0] > split_sizes[1], "We should have one biggest new split."
    return split_size_df["gpcsp"][0]


def read_gpcsp_probability_csv(in_path):
    "Read a CSV mapping GPCSPs to probabilities."
    return pd.read_csv(in_path, names=["gpcsp", "prob"])


def export_line_for_biggest_outside_gpcsp(original_path, outside_path, out_path):
    """Look for the outside GPCSP with the most taxa in the parent split, and export it
    and its probability to a CSV."""
    original_df = read_gpcsp_probability_csv(original_path)
    with_outside_df = read_gpcsp_probability_csv(outside_path)
    with_outside_gpcsps = set(with_outside_df["gpcsp"])
    original_gpcsps = set(original_df["gpcsp"])
    outside_gpcsps = list(with_outside_gpcsps.difference(original_gpcsps))
    if not outside_gpcsps:
        click.echo("No new GPCSPs found in " + outside_path)
        return
    to_export_df = with_outside_df.loc[
        with_outside_df["gpcsp"] == get_biggest_gpcsp(outside_gpcsps),
    ].copy()
    to_export_df["outside_split_count"] = len(outside_gpcsps)
    to_export_df.to_csv(out_path, index=False, header=False)
