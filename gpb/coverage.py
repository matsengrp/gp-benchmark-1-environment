"""functions for coverage analysis with mrbayes"""

import bito
import glob
import os
import time
import timeit
import numpy as np
import pandas as pd
import plotnine as p9
from csv import writer
from gpb.bitset_string import is_rootsplit
from gpb.ourbito import make_gp_instance

def posterior_credible_intervals(inst, gp_df, prefix):
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

    cols = ['pcsp', prefix + '_posterior', prefix + '_lb', prefix + '_ub', prefix + '_samples']

    posterior_df = pd.DataFrame(zip(pcsps, posterior_means, lbs, ubs, samples), 
                             columns = cols)

    return posterior_df


def run_estimation_benchmark(unif_newick_path, exp_newick_path, fasta_path, out_csv_prefix, tol, max_iter, use_gradients, mmap_path):

    inst = make_gp_instance(unif_newick_path, fasta_path, use_gradients, mmap_path)

    inst.estimate_branch_lengths(tol, max_iter, quiet = False, track_intermediate_iterations = False)

    inst.branch_lengths_to_csv(out_csv_prefix + ".gp_bl.csv")
    gp_df = pd.read_csv(out_csv_prefix + ".gp_bl.csv", header = None, names = ['pcsp', 'gp'])
    os.remove(out_csv_prefix + ".gp_bl.csv")
    unif_posterior_df = posterior_credible_intervals(inst, gp_df, prefix = 'unif')

    exp_inst = make_gp_instance(exp_newick_path, fasta_path, use_gradients, mmap_path)
    exp_posterior_df = posterior_credible_intervals(exp_inst, gp_df, prefix='exp')

    df = unif_posterior_df.merge(exp_posterior_df, on = 'pcsp').merge(gp_df, on = 'pcsp')
    df['dataset'] = out_csv_prefix

    df['gp_coverage'] = np.where(df['gp'] > df['unif_lb'],
                             np.where(df['gp'] < df['unif_ub'], 'covered', 'over'),
                             'under')
    df.to_csv(out_csv_prefix + ".csv")

    plot = (
        p9.ggplot(df, p9.aes(x = 'unif_posterior', y = 'gp'))
        + p9.geom_point()
        + p9.geom_abline(intercept = 0, slope = 1, color = 'blue')
    )
    plot.save(out_csv_prefix + ".pdf")


def run_timing_benchmark(newick_path, fasta_path, out_csv_prefix, tol, max_iter, use_gradients, mmap_path):
    start = time.time()

    inst = make_gp_instance(newick_path, fasta_path, use_gradients, mmap_path, print_stats = False)
    build_time_end = time.time()

    inst.estimate_branch_lengths(tol, max_iter, quiet = True, track_intermediate_iterations = False)
    estimation_time_end = time.time()

    os.remove(mmap_path)

    times = [build_time_end - start, estimation_time_end - build_time_end]

    outpath = out_csv_prefix + ".timing.csv"

    if not os.path.exists(outpath):
        headers = ['build_time', 'estimation_time']
        with open(outpath, 'w') as timingcsv:
            csvwriter = writer(timingcsv)
            csvwriter.writerow(headers)

    with open(outpath, "a+", newline = '') as timingcsv:
        csvwriter = writer(timingcsv)
        csvwriter.writerow(times)        


def run_coverage(datapath, sample_min):
    files = glob.glob(datapath + "/*.csv")
    dfs = []

    for filename in files:
        df = pd.read_csv(filename)
        dfs.append(df)

    full = pd.concat(dfs, axis = 0)
    full['dataset'] = full['dataset'].str.replace(".bench","")

    full['squared_error'] = (full['gp'] - full['unif_posterior'])**2
    mse = full.groupby(['dataset'])['squared_error'].mean().reset_index()

    r2 = full.groupby(['dataset'])[['unif_posterior', 'gp']].corr().iloc[0::2,-1].reset_index()[['dataset', 'gp']]
    r2['R-squared'] = r2['gp']**2

    summary_stats = r2.merge(mse, on = 'dataset')
    summary_stats.columns = ['dataset', 'corr', 'R-squared', 'MSE']
    summary_stats.to_csv(datapath + "/summary_stats.csv", index = False, float_format = '%.4f')

    greater_than = " ".join([">", str(sample_min), "samples"])
    less_than = " ".join(["<=", str(sample_min), "samples"])

    full['gtsamplemin'] = np.where(full['unif_samples'] >= sample_min, greater_than, less_than)

    coverage = full.groupby(['dataset', 'gtsamplemin', 'gp_coverage']).size().agg(
        {'count': lambda x: int(x), 'prop': lambda x: x / x.groupby(['dataset', 'gtsamplemin']).sum() * 100}
        ).unstack(level=0).reset_index()

    barplot = (
        p9.ggplot(coverage, p9.aes(x = 'dataset', y = 'prop', fill = 'gp_coverage', label = 'prop')) 
        + p9.geom_bar(stat=p9.stat_identity)
        + p9.geom_text(size = 5, position = p9.position_stack(vjust = 0.5),
                    format_string='{:.1f}%')
        + p9.facet_wrap('gtsamplemin')
        + p9.ylab('PCSP proportion')
        + p9.theme(axis_text_x = p9.element_text(size = 6))
    )
    barplot.save(datapath + "/coverage.pdf")


