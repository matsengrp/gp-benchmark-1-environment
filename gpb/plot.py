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

    df['iter'] = pd.to_numeric(df['iter'])
    df = df.loc[df.iter > 2,:]

    df['llh_diff'] = df.groupby('gpcsp')['per_pcsp_llh'].diff().fillna(0)
    gpcsps_to_check = df.loc[df.llh_diff < 0,:]['gpcsp'].unique()
    to_plot = df.loc[df['gpcsp'].isin(gpcsps_to_check).to_numpy().tolist()]

    plot = (
        p9.ggplot(to_plot,
                  p9.aes(x="iter", y = "llh_diff", group = 'gpcsp'))
        + p9.geom_point(p9.aes(color = 'gpcsp'))
        + p9.geom_path(p9.aes(color = 'gpcsp'))
        + p9.xlab("iterations")
        + p9.ylab("difference in per_pcsp_llh")
        + p9.ggtitle(out_path[0:3])
        + p9.theme(legend_position = 'none')
    )

    plot.save(filename = out_path)


def per_pcsp_likelihood_surfaces(surface_path, out_path):

    colnames = ['pcsp', 'branch_length', 'llh']
    df = pd.read_csv(surface_path, header = None, names = colnames)

    pcsps = pd.unique(df['pcsp'])

    plot_list = []

    for pcsp in pcsps:
        surface = df.loc[df.pcsp == pcsp]

        plot = (
            p9.ggplot(surface,
                      p9.aes(x='branch_length', y = 'llh'))
            + p9.geom_line()
            + p9.xlab('branch_length')
            + p9.ylab('per pcsp marginal log likelihood')
            + p9.ggtitle(pcsp)
            + p9.theme(plot_title = p9.element_text(size = 6))
        )
        
        plot_list.append(plot)

    p9.save_as_pdf_pages(plot_list, filename = out_path)


def per_pcsp_likelihood_surfaces_by_opt(nograd_surf_path, nograd_track_path, grad_surf_path, grad_track_path, out_path):

    colnames = ['gpcsp', 'branch_length', 'llh']
    
    nograd_surf = pd.read_csv(nograd_surf_path, header = None, names = colnames)
    nograd_track = pd.read_csv(nograd_track_path, header = None, names = colnames)
    nograd_surf['opt'] = 'brent'
    nograd_track['opt'] = 'brent'

    grad_surf = pd.read_csv(grad_surf_path, header = None, names = colnames)
    grad_track = pd.read_csv(grad_track_path, header = None, names = colnames)
    grad_surf['opt'] = 'newton'
    grad_track['opt'] = 'newton'

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
            + p9.geom_line()
            + p9.geom_point(p9.aes(x='branch_length', y='llh', color='iter'),
                            data = track)
            + p9.xlab('branch_length')
            + p9.ylab('per pcsp marginal log likelihood')
            + p9.facet_wrap('opt')
            + p9.ggtitle(gpcsp)
            + p9.theme(plot_title = p9.element_text(size = 6))
        )
        
        plot_list.append(plot)

    p9.save_as_pdf_pages(plot_list, filename = out_path)
