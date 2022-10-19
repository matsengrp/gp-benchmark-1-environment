"""functions for coverage analysis with mrbayes"""

import bito
import timeit
import glob
import os
import numpy as np
import pandas as pd
import plotnine as p9
from gpb.bitset_string import is_rootsplit
from gpb.ourbito import make_gp_instance

def posterior_credible_intervals(inst, gp_df, out_csv_prefix):
    """Get posterior credible intervals from MrBayes trees"""
    branch_lengths = {k: np.array(bl)
                      for k, bl in inst.gather_branch_lengths().items()
                      }

    sorted_keys = sorted(branch_lengths.keys())
    pcsps = []
    posterior_means = []
    lbs = []
    ubs = []
    samples = []

    for key in sorted_keys:
        pcsp = gp_df.loc[key]['pcsp']
        if is_rootsplit(pcsp):
            continue
        posterior_mean = np.mean(branch_lengths[key])
        lb = np.quantile(branch_lengths[key], 0.025)
        ub = np.quantile(branch_lengths[key], 0.975)
        num_samples = len(branch_lengths[key])
        pcsps.append(pcsp)
        posterior_means.append(posterior_mean)
        lbs.append(lb)
        ubs.append(ub)
        samples.append(num_samples)

    posterior_df = pd.DataFrame(zip(pcsps, posterior_means, lbs, ubs, samples), 
                             columns = ['pcsp', 'posterior', 'lb', 'ub', 'samples'])

    return posterior_df


def estimate_branch_lengths(inst, out_csv_prefix, tol, max_iter):
    inst.estimate_branch_lengths(tol, max_iter, quiet = False, track_intermediate_iterations = False)
    inst.branch_lengths_to_csv(out_csv_prefix + ".gp_bl.csv")
    gp_df = pd.read_csv(out_csv_prefix + ".gp_bl.csv", header = None, names = ['pcsp', 'gp'])
    os.remove(out_csv_prefix + ".gp_bl.csv")
    return gp_df


def run_benchmark(newick_path, fasta_path, out_csv_prefix, tol, max_iter, mmap_path):
    start = timeit.default_timer()
    use_gradients = True
    inst = make_gp_instance(newick_path, fasta_path, use_gradients, mmap_path)
    gp_df = estimate_branch_lengths(inst, out_csv_prefix, tol, max_iter)
    end = timeit.default_timer()
    time = end - start
    print("GP estimation took ", time, " seconds")

    # HN: This seems to be taken out of main. Not sure why or if it's needed.
    # inst.export_dag_to_json(out_csv_prefix + "_dag.json")
    
    posterior_df = posterior_credible_intervals(inst, gp_df, out_csv_prefix)

    df = posterior_df.merge(gp_df, on = 'pcsp')
    df['dataset'] = out_csv_prefix
    df['time'] = time

    df['coverage'] = np.where(df['gp'] > df['lb'],
                             np.where(df['gp'] < df['ub'], 'covered', 'over'),
                             'under')
    df.to_csv(out_csv_prefix + ".csv")

    plot = (
        p9.ggplot(df, p9.aes(x = 'posterior', y = 'gp'))
        + p9.geom_point()
        + p9.geom_abline(intercept = 0, slope = 1, color = 'blue')
    )
    plot.save(out_csv_prefix + ".pdf")


def run_coverage(datapath, sample_min):
    files = glob.glob(datapath + "/*.csv")
    dfs = []

    for filename in files:
        df = pd.read_csv(filename)
        dfs.append(df)

    full = pd.concat(dfs, axis = 0)
    full['dataset'] = full['dataset'].str.replace(".bench","")
    r2 = full.groupby(['dataset'])[['posterior', 'gp']].corr().iloc[0::2,-1].reset_index()[['dataset', 'gp']]
    r2['R-squared'] = r2['gp']**2
    r2.columns = ['dataset', 'corr', 'R-squared']
    r2.to_csv(datapath + "/r2.csv", index = False, float_format = '%.3f')

    greater_than = " ".join([">", str(sample_min), "samples"])
    less_than = " ".join(["<=", str(sample_min), "samples"])

    full['gtsamplemin'] = np.where(full['samples'] >= sample_min, greater_than, less_than)

    coverage = full.groupby(['dataset', 'gtsamplemin', 'coverage']).size().agg(
        {'count': lambda x: int(x), 'prop': lambda x: x / x.groupby(['dataset', 'gtsamplemin']).sum() * 100}
        ).unstack(level=0).reset_index()

    barplot = (
        p9.ggplot(coverage, p9.aes(x = 'dataset', y = 'prop', fill = 'coverage', label = 'prop')) 
        + p9.geom_bar(stat=p9.stat_identity)
        + p9.geom_text(size = 5, position = p9.position_stack(vjust = 0.5),
                    format_string='{:.1f}%')
        + p9.facet_wrap('gtsamplemin')
        + p9.ylab('PCSP proportion')
        + p9.theme(axis_text_x = p9.element_text(size = 6))
    )
    barplot.save(datapath + "/coverage.pdf")


