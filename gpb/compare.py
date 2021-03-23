"""Functions to compare parameter fits."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import logsumexp

import gpb.bitset_string as bitset_string

sns.set_context("poster")


def subsplit_probability_dict_of_csv(subsplit_csv):
    """Make a dictionary mapping from the subsplit (and its rotated version) to its
    probability given a subsplit probability CSV."""
    subsplit_df = pd.read_csv(subsplit_csv, names=["subsplit", "probability"])
    subsplit_dict = {
        row["subsplit"]: row["probability"] for _, row in subsplit_df.iterrows()
    }
    for subsplit, probability in subsplit_dict.copy().items():
        subsplit_dict[bitset_string.rotate_subsplit(subsplit)] = probability
    subsplit_dict[""] = 1.0
    return subsplit_dict


def compare_parameters(gp_csv, sa_csv, sa_subsplit_csv, out_prefix):
    """Compare parameters between GP and SA."""
    gp_name = "GP probability"
    sa_name = "SA probability"
    parent_probability_name = "parent probability"
    gp_df = pd.read_csv(gp_csv, names=["gpcsp", gp_name])
    sa_df = pd.read_csv(sa_csv, names=["gpcsp", sa_name])

    sa_probabilities = subsplit_probability_dict_of_csv(sa_subsplit_csv)
    sa_df["parent"] = sa_df["gpcsp"].apply(bitset_string.parent_bitset_of_gpcsp)
    sa_df["parent probability"] = sa_df["parent"].apply(sa_probabilities.get)

    df = pd.merge(gp_df, sa_df)
    # GP should have all of the SA GPCSPs, and also the fake ones.
    assert len(df) == len(sa_df)

    df["gpcsp"] = df["gpcsp"].apply(bitset_string.indexed_gpcsp_of_bitset_gpcsp)
    df["delta"] = df[sa_name] - df[gp_name]
    df = df[[parent_probability_name, "delta", gp_name, sa_name, "gpcsp"]]
    # Drop rows where GP and SA both have them as 1.
    df = df[(df[gp_name] < 1.0) | (df[sa_name] < 1.0)]
    df.to_csv(out_prefix + ".csv", index=False)

    correlation = df.corr().at[gp_name, sa_name]
    print(f"correlation: {correlation}")
    corr_df = pd.DataFrame({"dataset": [out_prefix], "correlation": [correlation]})
    corr_df.to_csv("correlation.csv", index=False)

    ax = sns.scatterplot(
        x=sa_name,
        y=gp_name,
        # If you want parent probability, add `size=parent_probability_name` here
        data=df,
        alpha=0.2,
    )
    ax.set_aspect(1)
    ax.set(ylim=ax.get_xlim())
    sns.despine()
    plt.legend([], [], frameon=False)
    plt.savefig(out_prefix + ".svg", bbox_inches="tight")


def compare_to_direct(direct_marginals_csv, prior_csv, sa_csv, out_prefix):
    """Compare SBN parameters estimated using the "direct" method to those from SA.

    The "direct" method here means the Bayes rule estimate using marginal likelihoods
    from just taking the mean from the trees (using GP branch lengths) found in the MB
    posterior sample.
    """
    direct_name = "Direct probability"
    sa_name = "SA probability"
    prior_df = pd.read_csv(prior_csv, names=["gpcsp", "prior"])
    df = pd.read_csv(direct_marginals_csv)
    gpcsp_count = len(df)
    df["gpcsp"] = df["gpcsp"].apply(bitset_string.pretty_gpcsp_of_string)
    df = pd.merge(df, prior_df)
    assert len(df) == gpcsp_count
    df["log_prior"] = df["prior"].apply(np.log)
    df["upost"] = df["marginal"] + df["log_prior"]
    df["parent"] = df.gpcsp.apply(bitset_string.parent_bitset_of_gpcsp)
    upost_sums = (
        df.groupby(["parent"])["upost"]
        .agg(logsumexp)
        .reset_index()
        .rename(columns={"upost": "upost_sum"})
    )
    df = pd.merge(upost_sums, df)
    assert len(df) == gpcsp_count
    df[direct_name] = (df["upost"] - df["upost_sum"]).apply(np.exp)
    sa_df = pd.read_csv(sa_csv, names=["gpcsp", sa_name])
    df = pd.merge(df, sa_df)
    assert len(df) == gpcsp_count
    df = df[(df[direct_name] < 1.0) | (df[sa_name] < 1.0)]
    correlation = df.corr().at[direct_name, sa_name]
    print(f"direct correlation: {correlation}")
    corr_df = pd.DataFrame(
        {"dataset": [out_prefix], "direct_correlation": [correlation]}
    )
    corr_df.to_csv("direct-correlation.csv", index=False)
    df.to_csv(out_prefix + ".csv", index=False)

    ax = sns.scatterplot(x=sa_name, y=direct_name, data=df, alpha=0.2)
    ax.set_aspect(1)
    ax.set(ylim=ax.get_xlim())
    sns.despine()
    plt.savefig(out_prefix + ".svg", bbox_inches="tight")


def add_metadata_to_sbn_df(df):
    """Add information about each GPCSP to a dataframe."""
    df["smaller_child_size"], df["larger_child_size"] = zip(
        *df["gpcsp"].apply(bitset_string.child_subsplit_taxon_set_sizes)
    )
    df["is_rootsplit"] = df["gpcsp"].apply(bitset_string.is_rootsplit)
    return df


def add_metadata_to_sbn_csv(sbn_csv_path, out_path):
    """Add information about each GPCSP to a CSV."""
    df = pd.read_csv(sbn_csv_path, names=["gpcsp", "prob"])
    add_metadata_to_sbn_df(df).to_csv(out_path, index=False)
