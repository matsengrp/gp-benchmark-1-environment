"""Functions to compare parameter fits."""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_context("poster")


def taxon_set_of_str_bitset(bitset):
    """Make the 1-indexed taxon set from a string bitset."""
    return {index + 1 for index, value in enumerate(bitset) if value == '1'}


def indexed_pcsp_of_bitset_pcsp(pcsp_str):
    """Make a string representation of a PCSP described in 1-indexed form."""
    [sister, focal, child1] = [taxon_set_of_str_bitset(bitset) for bitset in pcsp_str.split("|")]
    return str(((sister, focal), (set(sorted(focal - child1)), child1)))


def indexed_gpcsp_of_bitset_gpcsp(gpcsp_str):
    """Make a string representation of a GPCSP described in 1-indexed form."""
    if "|" in gpcsp_str:
        return indexed_pcsp_of_bitset_pcsp(gpcsp_str)
    # else:
    return str(taxon_set_of_str_bitset(gpcsp_str))


def flip_string_bits(bitset_string):
    """Exchange 0 and 1 in the bitset expressed as a string."""
    return bitset_string.translate(str.maketrans("01", "10"))


def sorted_parent_bitset_of_gpcsp(gpcsp_str):
    """Given a GPCSP string representation, make a sorted version of the parent
    subsplit."""
    if "|" in gpcsp_str:
        parent0, parent1, _ = gpcsp_str.split("|")
    else:
        parent0 = gpcsp_str
        parent1 = flip_string_bits(gpcsp_str)
    # TODO surprise! We have reverse sorting here!
    return "".join(sorted([parent0, parent1], reverse=True))


def sa_subsplit_probability_dict_of_csv(sa_subsplit_csv):
    """Make a dictionary mapping from the subsplit to the probability given a subsplit
    probability CSV."""
    sa_subsplit_df = pd.read_csv(sa_subsplit_csv, names=["subsplit", "probability"], header=None)
    sa_subsplit_dict = {row["subsplit"]: row["probability"] for _, row in sa_subsplit_df.iterrows()}
    sa_subsplit_dict[""] = 1.
    return sa_subsplit_dict


sa_subsplit_csv = "/loc/matsen/gp-benchmark-1/ds1/_output/sa-sbn-parameters.subsplit.csv"
sa_csv = "/loc/matsen/gp-benchmark-1/ds1/_output/sa-sbn-parameters.csv"
sa_name = "SA probability"
sa_df = pd.read_csv(sa_csv, names=["gpcsp", sa_name])

sa_probabilities = sa_subsplit_probability_dict_of_csv(sa_subsplit_csv)
sa_df["sorted parent"] = sa_df["gpcsp"].apply(sorted_parent_bitset_of_gpcsp)
sa_df["parent_probability"] = sa_df["sorted parent"].apply(sa_probabilities.get)
sa_df.head()

# TODO what's up with the nans here?

for _, row in sa_df.loc[sa_df["parent_probability"].isnull()].iterrows():
    print(row)


gpcsp = sa_df.iloc[4,0]

pd.read_csv(sa_subsplit_csv)


def compare_parameters(gp_csv, sa_csv, out_prefix):
    """Compare parameters between GP and SA."""
    gp_name = "GP probability"
    sa_name = "SA probability"
    # TODO note that we are using the rootsplit as a fake header row?
    gp_df = pd.read_csv(gp_csv, names=["gpcsp", gp_name])
    sa_df = pd.read_csv(sa_csv, names=["gpcsp", sa_name])
    df = pd.merge(gp_df, sa_df)
    # GP should have all of the SA GPCSPs, and also the fake ones.
    assert len(df) == len(sa_df)

    df["gpcsp"] = df["gpcsp"].apply(indexed_gpcsp_of_bitset_gpcsp)
    df["delta"] = df[sa_name] - df[gp_name]
    df = df[["delta", gp_name, sa_name, "gpcsp"]]
    # Drop rows where GP and SA both have them as 1.
    df = df[(df[gp_name] < 1.) | (df[sa_name] < 1.)]
    df.to_csv(out_prefix + ".csv", index=False)

    correlation = df.corr().at[gp_name, sa_name]
    print(f"correlation: {correlation}")
    corr_df = pd.DataFrame({"dataset": [out_prefix], "correlation": [correlation]})
    corr_df.to_csv("correlation.csv", index=False)

    ax = sns.scatterplot(x=sa_name, y=gp_name, data=df, alpha=0.2)
    ax.set_aspect(1)
    ax.set(ylim=ax.get_xlim())
    sns.despine()
    plt.savefig(out_prefix + ".svg", bbox_inches="tight")
