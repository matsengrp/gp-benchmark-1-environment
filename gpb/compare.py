"""
Functions to compare parameter fits.
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_context("poster")


def compare_parameters(gp_csv, sa_csv, out_prefix):
    """
    Compare parameters between GP and SA.
    """
    gp_df = pd.read_csv(gp_csv, names=["gpcsp", "GP probability"])
    sa_df = pd.read_csv(sa_csv, names=["gpcsp", "SA probability"])
    df = pd.merge(gp_df, sa_df)
    # GP should have all of the SA GPCSPs, and also the fake ones.
    assert len(df) == len(sa_df)
    df.to_csv(out_prefix + ".csv", index=False)

    print(df.corr())

    ax = sns.scatterplot(x="SA probability", y="GP probability", data=df, alpha=0.2)
    ax.set_aspect(1)
    ax.set(ylim=ax.get_xlim())
    sns.despine()
    plt.savefig(out_prefix + ".svg", bbox_inches="tight")
