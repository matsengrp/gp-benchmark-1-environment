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


def parent_bitset_of_gpcsp(gpcsp_str):
    """Given a PCSP string representation, extract the parent subsplit with no bar.
    Return empty string if there is no bar in the input (a rootsplit)."""
    if "|" in gpcsp_str:
        parent0, parent1, _ = gpcsp_str.split("|")
        return "".join([parent0, parent1])
    # else:
    return ""


def rotate_subsplit(bitset_string):
    """Given a bitset abcdef, return defabc."""
    assert len(bitset_string) % 2 == 0
    half_len = len(bitset_string) // 2
    return bitset_string[half_len:] + bitset_string[:half_len]


def sa_subsplit_probability_dict_of_csv(sa_subsplit_csv):
    """Make a dictionary mapping from the subsplit (and its rotated version) to its
    probability given a subsplit probability CSV."""
    sa_subsplit_df = pd.read_csv(sa_subsplit_csv, names=["subsplit", "probability"])
    sa_subsplit_dict = {row["subsplit"]: row["probability"] for _, row in sa_subsplit_df.iterrows()}
    for subsplit, probability in sa_subsplit_dict.copy().items():
        sa_subsplit_dict[rotate_subsplit(subsplit)] = probability
    sa_subsplit_dict[""] = 1.
    return sa_subsplit_dict


def compare_parameters(gp_csv, sa_csv, sa_subsplit_csv, out_prefix):
    """Compare parameters between GP and SA."""
    gp_name = "GP probability"
    sa_name = "SA probability"
    parent_probability_name = "parent probability"
    gp_df = pd.read_csv(gp_csv, names=["gpcsp", gp_name])
    sa_df = pd.read_csv(sa_csv, names=["gpcsp", sa_name])

    sa_probabilities = sa_subsplit_probability_dict_of_csv(sa_subsplit_csv)
    sa_df["parent"] = sa_df["gpcsp"].apply(parent_bitset_of_gpcsp)
    sa_df["parent probability"] = sa_df["parent"].apply(sa_probabilities.get)

    df = pd.merge(gp_df, sa_df)
    # GP should have all of the SA GPCSPs, and also the fake ones.
    assert len(df) == len(sa_df)

    df["gpcsp"] = df["gpcsp"].apply(indexed_gpcsp_of_bitset_gpcsp)
    df["delta"] = df[sa_name] - df[gp_name]
    df = df[[parent_probability_name, "delta", gp_name, sa_name, "gpcsp"]]
    # Drop rows where GP and SA both have them as 1.
    df = df[(df[gp_name] < 1.) | (df[sa_name] < 1.)]
    df.to_csv(out_prefix + ".csv", index=False)

    correlation = df.corr().at[gp_name, sa_name]
    print(f"correlation: {correlation}")
    corr_df = pd.DataFrame({"dataset": [out_prefix], "correlation": [correlation]})
    corr_df.to_csv("correlation.csv", index=False)

    ax = sns.scatterplot(x=sa_name, y=gp_name, size=parent_probability_name, data=df, alpha=0.2)
    ax.set_aspect(1)
    ax.set(ylim=ax.get_xlim())
    sns.despine()
    plt.legend([],[], frameon=False)
    plt.savefig(out_prefix + ".svg", bbox_inches="tight")
