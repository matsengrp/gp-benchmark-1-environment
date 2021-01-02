"""Functions to compare parameter fits."""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from scipy.special import logsumexp

sns.set_context("poster")


def taxon_set_of_str_bitset(bitset):
    """Make the 1-indexed taxon set from a string bitset."""
    return {index + 1 for index, value in enumerate(bitset) if value == "1"}


def indexed_pcsp_of_bitset_pcsp(pcsp_str):
    """Make a string representation of a PCSP described in 1-indexed form."""
    [sister, focal, child1] = [
        taxon_set_of_str_bitset(bitset) for bitset in pcsp_str.split("|")
    ]
    return str(((sister, focal), (set(sorted(focal - child1)), child1)))


def indexed_gpcsp_of_bitset_gpcsp(gpcsp_str):
    """Make a string representation of a pretty GPCSP described in 1-indexed form."""
    if "|" in gpcsp_str:
        return indexed_pcsp_of_bitset_pcsp(gpcsp_str)
    # else:
    # We have a rootsplit.
    return str(taxon_set_of_str_bitset(gpcsp_str))


def parent_bitset_of_gpcsp(gpcsp_str, sort=False):
    """Given a PCSP string representation, extract the parent subsplit with no bar.
    Return empty string if there is no bar in the input (a rootsplit)."""
    if "|" in gpcsp_str:
        parent0, parent1, _ = gpcsp_str.split("|")
        subsplits = [parent0, parent1]
        if sort:
            subsplits.sort()
        return "".join(subsplits)
    # else:
    return ""


def rotate_subsplit(bitset_string):
    """Given a bitset abcdef, return defabc."""
    assert len(bitset_string) % 2 == 0
    half_len = len(bitset_string) // 2
    return bitset_string[half_len:] + bitset_string[:half_len]


def pretty_gpcsp_of_string(s):
    """Insert vertical bars to make a PCSP bitset string easier to read."""
    length = len(s)
    assert length % 3 == 0
    chunk_length = length // 3
    return "|".join([s[i : i + chunk_length] for i in range(0, length, chunk_length)])


def subsplit_probability_dict_of_csv(subsplit_csv):
    """Make a dictionary mapping from the subsplit (and its rotated version) to its
    probability given a subsplit probability CSV."""
    subsplit_df = pd.read_csv(subsplit_csv, names=["subsplit", "probability"])
    subsplit_dict = {
        row["subsplit"]: row["probability"] for _, row in subsplit_df.iterrows()
    }
    for subsplit, probability in subsplit_dict.copy().items():
        subsplit_dict[rotate_subsplit(subsplit)] = probability
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
    sa_df["parent"] = sa_df["gpcsp"].apply(parent_bitset_of_gpcsp)
    sa_df["parent probability"] = sa_df["parent"].apply(sa_probabilities.get)

    df = pd.merge(gp_df, sa_df)
    # GP should have all of the SA GPCSPs, and also the fake ones.
    assert len(df) == len(sa_df)

    df["gpcsp"] = df["gpcsp"].apply(indexed_gpcsp_of_bitset_gpcsp)
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
    df["gpcsp"] = df["gpcsp"].apply(pretty_gpcsp_of_string)
    df = pd.merge(df, prior_df)
    assert len(df) == gpcsp_count
    df["log_prior"] = df["prior"].apply(np.log)
    df["upost"] = df["marginal"] + df["log_prior"]
    df["parent"] = df.gpcsp.apply(parent_bitset_of_gpcsp)
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
    print(df.corr().at[direct_name, sa_name])
    df.to_csv(out_prefix + ".csv", index=False)

    ax = sns.scatterplot(x=sa_name, y=direct_name, data=df, alpha=0.2)
    ax.set_aspect(1)
    ax.set(ylim=ax.get_xlim())
    sns.despine()
    plt.savefig(out_prefix + ".svg", bbox_inches="tight")
