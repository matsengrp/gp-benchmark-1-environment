"""Functions to compare parameter fits."""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_context("poster")


def compare_parameters(gp_csv, sa_csv, out_prefix):
    """Compare parameters between GP and SA."""
    gp_name = "GP probability"
    sa_name = "SA probability"
    gp_df = pd.read_csv(gp_csv, names=["gpcsp", gp_name])
    sa_df = pd.read_csv(sa_csv, names=["gpcsp", sa_name])
    df = pd.merge(gp_df, sa_df)
    # GP should have all of the SA GPCSPs, and also the fake ones.
    assert len(df) == len(sa_df)
    df.to_csv(out_prefix + ".csv", index=False)

    correlation = df.corr().at[gp_name, sa_name]
    print("correlation:" + correlation)
    corr_df = pd.DataFrame({"dataset": [out_prefix], "correlation": [correlation]})
    corr_df.to_csv("correlation.csv", index=False)

    ax = sns.scatterplot(x=sa_name, y=gp_name, data=df, alpha=0.2)
    ax.set_aspect(1)
    ax.set(ylim=ax.get_xlim())
    sns.despine()
    plt.savefig(out_prefix + ".svg", bbox_inches="tight")
