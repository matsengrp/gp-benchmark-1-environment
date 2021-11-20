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

def per_pcsp_likelihood_surfaces_by_opt(nograd_surf_path, nograd_track_path, grad_surf_path, grad_track_path, out_path):

    colnames = ['gpcsp', 'branch_length', 'llh']
    
    nograd_surf = pd.read_csv(nograd_surf_path, header = None, names = colnames)
    nograd_track = pd.read_csv(nograd_track_path, header = None, names = colnames)
    nograd_surf['opt'] = 'nograd'
    nograd_track['opt'] = 'nograd'

    grad_surf = pd.read_csv(grad_surf_path, header = None, names = colnames)
    grad_track = pd.read_csv(grad_track_path, header = None, names = colnames)
    grad_surf['opt'] = 'grad'
    grad_track['opt'] = 'grad'

    surf_df = nograd_surf.append(grad_surf, ignore_index=True, sort=False)
    track_df = nograd_track.append(grad_track, ignore_index=True, sort=False)

    track_df = track_df.drop_duplicates()
    track_df['iter'] = pd.Categorical((track_df.groupby(['opt','gpcsp']).cumcount()) +1)

    gpcsps = pd.unique(surf_df['gpcsp'])

    plot_list = []

    for gpcsp in gpcsps:
        surface = surf_df[surf_df['gpcsp'] == gpcsp]
        track = track_df[track_df['gpcsp'] == gpcsp]

        plot = (
            p9.ggplot(surface,
                      p9.aes(x='branch_length', y = 'llh'))
            + p9.geom_line(p9.aes(linetype= 'opt'))
            + p9.geom_point(p9.aes(x='branch_length', y='llh',color='iter'),
                            data = track)
            + p9.xlab('branch_length')
            + p9.ylab('per pcsp marginal log likelihood')
            + p9.ggtitle(gpcsp)
            + p9.theme(plot_title = p9.element_text(size = 6))
        )
        
        plot_list.append(plot)

    p9.save_as_pdf_pages(plot_list, filename = out_path)
