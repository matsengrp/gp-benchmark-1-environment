"""Functions to compare subsplits outside the support."""

import click
import pandas as pd

import gpb.bitset_string as bitset_string
import gpb.compare as compare


def get_biggest_gpcsp(gpcsp_list):
    """Return the GPCSP from the list with the most taxa in its parent split."""
    split_size_df = (
        pd.DataFrame(
            {
                "gpcsp": gpcsp_list,
                "parent_split_size": [
                    bitset_string.parent_split_size(gpcsp) for gpcsp in gpcsp_list
                ],
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
    and information about it to a CSV."""
    original_df = read_gpcsp_probability_csv(original_path)
    with_outside_df = read_gpcsp_probability_csv(outside_path)
    with_outside_gpcsps = set(with_outside_df["gpcsp"])
    original_gpcsps = set(original_df["gpcsp"])
    outside_gpcsps = list(with_outside_gpcsps.difference(original_gpcsps))
    if not outside_gpcsps:
        click.echo(
            f"No new GPCSPs found in {outside_path}; this doesn't indicate an error "
            "but explains why aren't making an `.outside.csv` in this case."
        )
        return
    to_export_df = with_outside_df.loc[
        with_outside_df["gpcsp"] == get_biggest_gpcsp(outside_gpcsps),
    ].copy()
    to_export_df["outside_split_count"] = len(outside_gpcsps)
    compare.add_metadata_to_sbn_df(to_export_df).to_csv(out_path, index=False)
