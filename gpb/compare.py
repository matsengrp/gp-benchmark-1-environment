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


def compare_parameters(gp_csv, sa_csv, out_prefix):
    """Compare parameters between GP and SA."""
    gp_name = "GP probability"
    sa_name = "SA probability"
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
