"""Functions to compare subsplits outside the support."""

import click
import pandas as pd


def split_pcsp_into_list(pcsp):
    """Split a "pretty" string representation of a PCSP into three components."""
    pcsp_list = pcsp.split("|")
    assert len(pcsp_list) == 3
    return pcsp_list


def parent_split_size(gpcsp):
    """Count the number of taxa in the parent split of a PCSP, or return 1+taxon_count if
    we have a rootsplit."""
    if "|" in gpcsp:
        pcsp_list = split_pcsp_into_list(gpcsp)
        return (pcsp_list[0] + pcsp_list[1]).count("1")
    # else we have a rootsplit
    return 1 + len(gpcsp)


def child_subsplit_taxon_set_sizes(gpcsp):
    """Given a "pretty" string representation of a PCSP, find the number of taxa in the
    smallest component of the child subsplit."""
    if "|" in gpcsp:
        pcsp_list = split_pcsp_into_list(gpcsp)
        child_component_size_2 = pcsp_list[2].count("1")
        child_component_size_1 = pcsp_list[1].count("1") - child_component_size_2
    else:
        child_component_size_1 = gpcsp.count("0")
        child_component_size_2 = gpcsp.count("1")
    return (
        min(child_component_size_1, child_component_size_2),
        max(child_component_size_1, child_component_size_2),
    )


def get_biggest_gpcsp(gpcsp_list):
    """Return the GPCSP from the list with the most taxa in its parent split."""
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
    chosen_gpcsp = to_export_df.iloc[0, 0]
    to_export_df["outside_split_count"] = len(outside_gpcsps)
    smaller_child_size, larger_child_size = child_subsplit_taxon_set_sizes(chosen_gpcsp)
    to_export_df["smaller_child_size"] = smaller_child_size
    to_export_df["larger_child_size"] = larger_child_size
    to_export_df["is_rootsplit"] = "|" not in chosen_gpcsp
    to_export_df.to_csv(out_path, index=False)
