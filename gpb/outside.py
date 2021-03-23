"""Functions to compare subsplits outside the support."""

import click
import pandas as pd
import plotnine as p9

import gpb.bitset_string as bitset_string
import gpb.compare as compare

p9.theme_set(p9.theme_bw())


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


def plot_outside_prob_comparison(outside_csv_path, inside_csv_path, out_path):
    """Make a box and whisker plot comparing probabilities inside and outside of the support."""
    df = pd.read_csv(outside_csv_path)
    df.drop("gpcsp", axis=1, inplace=True)
    df["in_support"] = False

    inside_df = pd.read_csv(inside_csv_path)
    inside_df.drop("gpcsp", axis=1, inplace=True)
    inside_df = inside_df[inside_df["prob"] < 1.0]
    inside_df["in_support"] = True

    df = df.merge(inside_df, how="outer").reset_index()
    df["is_rootsplit"].sum()
    # We are fixing the outgroup, so aren't interested in NNIs that modify the rootsplit.
    df = df[~df["is_rootsplit"].astype(bool)]
    df.drop("is_rootsplit", axis=1, inplace=True)

    plot = (
        p9.ggplot(df)
        + p9.geom_boxplot(
            p9.aes(x="factor(smaller_child_size)", y="prob", color="in_support")
        )
        + p9.scale_y_log10()
        + p9.scale_color_brewer(type="qual", palette="Dark2")
        + p9.xlab("size of the smaller child clade")
        + p9.ylab("normalized probability via composite likelihood")
    )
    plot.save(out_path)
