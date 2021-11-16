"""Plotting per pcsp outputs from optimization"""

import pandas as pd
import plotnine as p9

import gpb.bitset_string as bitset_string
import gpb.compare as compare

p9.theme_set(p9.theme_bw())

def per_pcsp_likelihoods_from_opt_plot(per_pcsp_likelihoods_path, out_path):

    df = pd.read_csv(per_pcsp_likelihoods_path, header = None)
    df = df.rename(columns = {0:'gpcsp'})
    
    df = compare.add_metadata_to_sbn_df(df)
    df = df.drop(['larger_child_size'], axis = 1)

    df = pd.melt(df, id_vars = ['gpcsp', 'smaller_child_size', 'is_rootsplit'], 
                 var_name='iter', value_name='per_pcsp_llh')

    plot = (
        p9.ggplot(df,
                  p9.aes(x="iter", y = "per_pcsp_llh", group = "gpcsp"))
        + p9.geom_line(p9.aes(color ="factor(smaller_child_size)"))   
        + p9.xlab("iterations over composite marginal convergence")
        + p9.ylab("per pcsp marginal log likelihood")
    )
    plot.save(out_path)

def per_pcsp_likelihood_surfaces(per_pcsp_likelihood_surfaces_path, out_path):

    df = pd.read_csv(per_pcsp_likelihood_surfaces_path, header = None)
    df = df.rename(columns = {0:'gpcsp', 1:'branch_length', 2:'llh'})

    df = df[df['gpcsp'] != df['gpcsp'][0]]

    df = compare.add_metadata_to_sbn_df(df)
    
    df = df.drop(['larger_child_size'], axis = 1)

    df['llh'] = df.groupby('gpcsp')['llh'].transform(
        lambda x: 0 if (x.min() == x.max()) else ((x - x.min())/(x.max()-x.min()))
    )

    plot = (
        p9.ggplot(df,
                  p9.aes(x="branch_length", y = "llh", group = "gpcsp"))
        + p9.geom_line()
        + p9.facets.facet_wrap('smaller_child_size', labeller = 'label_both') 
        + p9.xlab("branch length")
        + p9.ylab("per pcsp marginal log likelihood")
    )
    plot.save(out_path)
