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
        try: 
            pcsp = gp_df.loc[key]['pcsp']
        except:
            continue

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
    exp_inst.estimate_branch_lengths(tol, max_iter, quiet = True, track_intermediate_iterations = False)
    exp_inst.branch_lengths_to_csv(out_csv_prefix + ".exp_gp_bl.csv")
    exp_gp_df = pd.read_csv(out_csv_prefix + ".exp_gp_bl.csv", header = None, names = ['pcsp', 'gp'])
    os.remove(out_csv_prefix + ".exp_gp_bl.csv")
    exp_posterior_df = posterior_credible_intervals(exp_inst, exp_gp_df, prefix='exp')

    df = unif_posterior_df.merge(exp_posterior_df, how = 'outer', on = 'pcsp').merge(gp_df, how = 'left', on = 'pcsp')
    df['dataset'] = out_csv_prefix

    df['gp_coverage'] = np.where(df['gp'] > df['unif_lb'],
                             np.where(df['gp'] < df['unif_ub'], 'covered', 'over'),
                             'under')
    df.to_csv(out_csv_prefix + ".csv")


def run_timing_benchmark(newick_path, fasta_path, out_csv_prefix, tol, max_iter, use_gradients, mmap_path):
    start = time.time()

    inst = make_gp_instance(newick_path, fasta_path, use_gradients, mmap_path, print_stats = False)
    build_time_end = time.time()

    inst.estimate_branch_lengths(tol, max_iter, quiet = True, track_intermediate_iterations = False)
    estimation_time_end = time.time()

    os.remove(mmap_path)

    row = [out_csv_prefix, build_time_end - start, estimation_time_end - build_time_end]

    outpath = out_csv_prefix + ".csv"

    if not os.path.exists(outpath):
        headers = ['dataset', 'build_time', 'estimation_time']
        with open(outpath, 'w') as timingcsv:
            csvwriter = writer(timingcsv)
            csvwriter.writerow(headers)

    with open(outpath, "a+", newline = '') as timingcsv:
        csvwriter = writer(timingcsv)
        csvwriter.writerow(row)        


def compile_timing_stats(datapath, uniq = False):
    if uniq:
        suff = 'uniq'
        files = glob.glob(datapath + "/*.timing.uniq.csv")
    else:
        suff = 'full'
        files = glob.glob(datapath + "/*.timing.csv")

    dfs = []
    for filename in files:
        df = pd.read_csv(filename)
        dfs.append(df)

    full = pd.concat(dfs, axis = 0)
    full['dataset'] = full['dataset'].str.replace(r'.bench(.*)','', regex = True)

    time_stats = full.groupby(['dataset'])[['build_time', 'estimation_time']].aggregate(['mean', 'std']).reset_index()

    cols = ['dataset', 'build_mean_' + suff, 'build_std_' + suff, 'estimation_mean_' + suff, 'estimation_std_' + suff]

    time_stats.columns = cols
    return(time_stats)


def convert_pcsp(pcsp_string, input = 'vbpi', output = 'gp', ds = 'ds1'):
    gp_taxon_order = {
            'ds1': [1, 24, 10, 22, 25, 7, 11, 16, 20, 18, 12, 4, 17, 6, 8, 13, 9, 3, 14, 19, 21, 5, 2, 23, 26, 27, 15], 
            'ds2': [], 
            'ds3': [1, 31, 14, 2, 34, 8, 10, 35, 4, 12, 15, 30, 32, 19, 3, 11, 18, 13, 16, 9, 36, 33, 20, 23, 25, 26, 27, 28, 24, 21, 22, 29, 5, 17, 6, 7], 
            'ds4': [1, 4, 11, 10, 37, 22, 9, 14, 36, 41, 17, 31, 32, 15, 16, 39, 33, 6, 2, 23, 13, 7, 18, 27, 28, 19, 20, 3, 35, 21, 29, 26, 40, 30, 38, 12, 24, 34, 5, 25, 8], 
            'ds5': [1, 10, 39, 36, 5, 27, 45, 35, 24, 13, 34, 47, 48, 12, 33, 38, 30, 4, 43, 14, 49, 46, 20, 42, 22, 28, 32, 11, 15, 50, 8, 40, 7, 25, 17, 6, 23, 9, 31, 16, 21, 26, 29, 37, 41, 44, 18, 19, 2, 3], 
            'ds6': [1, 46, 10, 38, 47, 9, 11, 28, 29, 37, 12, 25, 36, 41, 7, 8, 13, 14, 3, 22, 23, 31, 48, 16, 42, 17, 43, 19, 26, 27, 18, 32, 44, 4, 2, 20, 5, 35, 21, 24, 34, 33, 45, 39, 40, 15, 30, 49, 50, 6], 
            'ds7': [1, 53, 20, 2, 56, 8, 16, 4, 57, 18, 21, 52, 54, 25, 3, 10, 11, 9, 14, 15, 12, 13, 19, 22, 58, 59, 26, 29, 27, 28, 32, 33, 43, 44, 30, 31, 47, 48, 49, 50, 45, 46, 34, 35, 36, 37, 38, 39, 40, 42, 41, 51, 5, 55, 23, 17, 24, 6, 7], 
            'ds8': [1, 17, 18, 19, 32, 33, 34, 51, 37, 7, 15, 60, 61, 11, 12, 16, 50, 3, 27, 29, 58, 59, 30, 31, 28, 13, 8, 35, 39, 40, 44, 45, 63, 64, 55, 62, 9, 38, 54, 14, 2, 36, 20, 25, 26, 56, 57, 21, 23, 22, 24, 47, 48, 49, 52, 53, 42, 43, 46, 10, 6, 41, 5, 4]
    }
    taxon_count = {ds: len(gp_taxon_order[ds]) for ds in gp_taxon_order}
    rimu_taxon_order = {ds: list(range(1, taxon_count[ds] + 1)) for ds in taxon_count}        
    vbpi_taxon_order = {
            'ds1': [15, 18, 20, 16, 11, 25, 7, 22, 10, 24, 1, 27, 12, 4, 8, 17, 6, 21, 19, 14, 13, 9, 3, 5, 26, 23, 2],
            'ds2': [],
            'ds3': [7, 6, 18, 11, 17, 33, 36, 16, 13, 9, 29, 22, 21, 24, 28, 27, 26, 25, 23, 20, 5, 19, 32, 30, 15, 12, 3, 10, 35, 4, 34, 8, 2, 31, 1, 14],
            'ds4': [8, 24, 12, 34, 5, 25, 38, 30, 11, 4, 1, 39, 16, 15, 40, 28, 27, 18, 26, 29, 21, 20, 19, 35, 3, 37, 22, 10, 9, 33, 6, 32, 31, 17, 41, 36, 14, 23, 2, 13, 7],
            'ds5': [1, 3, 2, 19, 18, 44, 41, 29, 37, 26, 21, 16, 31, 11, 15, 25, 50, 8, 40, 7, 46, 28, 22, 42, 20, 49, 14, 17, 6, 32, 38, 33, 12, 48, 47, 34, 13, 24, 35, 45, 27, 36, 39, 10, 5, 30, 43, 4, 9, 23],
            'ds6': [6, 46, 1, 37, 29, 28, 11, 25, 12, 41, 36, 8, 7, 50, 49, 30, 15, 10, 47, 38, 9, 40, 39, 45, 33, 34, 24, 42, 16, 48, 31, 23, 13, 22, 14, 3, 18, 27, 26, 43, 19, 17, 44, 32, 4, 20, 2, 35, 21, 5],
            'ds7': [7, 6, 25, 54, 52, 21, 18, 3, 16, 57, 4, 56, 8, 2, 20, 53, 1, 24, 17, 23, 55, 59, 58, 22, 19, 14, 15, 13, 12, 10, 11, 9, 5, 51, 46, 45, 50, 49, 48, 47, 31, 30, 44, 43, 33, 32, 28, 27, 29, 26, 37, 36, 35, 34, 41, 38, 42, 40, 39],
            'ds8': [4, 5, 41, 10, 6, 46, 43, 42, 36, 28, 31, 30, 59, 58, 29, 27, 50, 16, 12, 11, 3, 38, 62, 55, 9, 64, 63, 45, 44, 40, 39, 35, 13, 8, 54, 14, 2, 61, 60, 15, 51, 34, 33, 32, 37, 19, 7, 18, 17, 1, 53, 52, 49, 48, 47, 24, 22, 23, 21, 57, 56, 26, 25, 20]
    }

    gp_to_rimu = {ds: [gp_taxon_order[ds].index(taxon) for taxon in rimu_taxon_order[ds]] for ds in gp_taxon_order}
    gp_to_vbpi = {ds: [gp_taxon_order[ds].index(taxon) for taxon in vbpi_taxon_order[ds]] for ds in gp_taxon_order}
    rimu_to_gp = {ds: [rimu_taxon_order[ds].index(taxon) for taxon in gp_taxon_order[ds]] for ds in rimu_taxon_order}
    rimu_to_vbpi = {ds: [rimu_taxon_order[ds].index(taxon) for taxon in vbpi_taxon_order[ds]] for ds in rimu_taxon_order}
    vbpi_to_gp = {ds: [vbpi_taxon_order[ds].index(taxon) for taxon in gp_taxon_order[ds]] for ds in vbpi_taxon_order}
    vbpi_to_rimu = {ds: [vbpi_taxon_order[ds].index(taxon) for taxon in rimu_taxon_order[ds]] for ds in vbpi_taxon_order}
    type_to_type = {
            "gp": {"rimu": gp_to_rimu, "vbpi": gp_to_vbpi},
            "rimu": {"gp": rimu_to_gp, "vbpi": rimu_to_vbpi},
            "vbpi": {"gp": vbpi_to_gp, "rimu": vbpi_to_rimu}
    }
    input_to_output = type_to_type[input][output][ds]
    mid_point = taxon_count[ds]
    clade_reorder = lambda bits: "".join([bits[input_to_output[j]] for j in range(mid_point)])
    choice = lambda p, c: min(c, bin(int(p, 2) & ~int(c, 2))[2:].zfill(mid_point))

    if len(pcsp_string) == mid_point:
        return clade_reorder(pcsp_string)
    else:
        sister = clade_reorder(pcsp_string[:mid_point])
        focal = clade_reorder(pcsp_string[mid_point : 2 * mid_point])
        child = clade_reorder(pcsp_string[2 * mid_point :])
        return sister + focal + choice(focal, child)


def merge_gp_vbpi(datapath, ds):
    gp_path = glob.glob(datapath + '/*' + ds + '*.estimation.csv')[0]
    vbpi_path = glob.glob(datapath + '/vbpi/*' + ds + '*.csv')[0]

    gp = pd.read_csv(gp_path)
    gp['pcsp_nobar'] = gp.pcsp.str.replace("|", "", regex = False)

    vbpi_cols = ['vbpi_pcsp', 'vbpi']
    vbpi = pd.read_csv(vbpi_path, names = vbpi_cols)
    vbpi['vbpi_pcsp_nobar'] = vbpi.vbpi_pcsp.str.replace("|", "", regex = False)
    vbpi['pcsp_nobar'] = vbpi.vbpi_pcsp_nobar.apply(lambda x: convert_pcsp(x, input = 'vbpi', output = 'gp', ds = ds))

    merged = gp.merge(vbpi, how = 'outer', on = 'pcsp_nobar')
    merged['dataset'] = ds

    merged['vbpi_coverage'] = np.where(merged['vbpi'] > merged['exp_lb'],
                                       np.where(merged['vbpi'] < merged['exp_ub'], 
                                                'covered', 'over'),'under')
    merged.to_csv(datapath + '/' + ds + '_merged.csv')
    return merged


def get_estimation_stats(df, estimate, truth, coverage, out_prefix):
    label = df['dataset'].iloc[0]
    corr = df[[truth, estimate]].corr().iloc[0::2, -1].item()
    r2 = corr**2
    mae = (abs(df[truth] - df[estimate])).mean()
    coverage_stat = np.where(df[coverage] == 'covered', 1, 0).mean()
    output = {
        'dataset': label,
        estimate+'_correlation': corr,
        estimate+'_r_squared': r2,
        estimate+'_mae': mae,
        estimate+'_coverage': coverage_stat
    }

    plot = (
        p9.ggplot(df, p9.aes(x = truth, y = estimate))
        + p9.geom_point()
        + p9.geom_abline(intercept = 0, slope = 1, color = 'blue')
        + p9.ggtitle(label + ': Correlation = ' + str(corr))
    )
    plot.save(out_prefix + ".pdf")
    return output


def compile_estimation_stats(datapath, sample_min):
    ds_list = ['ds1', 'ds3', 'ds4', 'ds5', 'ds6', 'ds7', 'ds8']

    gp_stats_full = []
    vbpi_stats_full  = []

    gp_stats_samplesubset = []
    vbpi_stats_samplesubset = []

    for ds in ds_list:
        df = merge_gp_vbpi(datapath, ds)

        gp_df = df.loc[(~df['unif_samples'].isna()) & (~df['gp'].isna())]
        gp_fulldf_prefix = ds + '_gp_full'
        gp_stats_full.append(get_estimation_stats(gp_df, 'gp', 'unif_posterior', 'gp_coverage', gp_fulldf_prefix))

        vbpi_df = df.loc[(~df['unif_samples'].isna()) & (~df['vbpi'].isna())]
        vbpi_fulldf_prefix = ds + '_vbpi_full'
        vbpi_stats_full.append(get_estimation_stats(vbpi_df,'vbpi', 'exp_posterior', 'vbpi_coverage', vbpi_fulldf_prefix))

        gp_subset = gp_df.loc[gp_df['unif_samples'] >= sample_min]
        gp_subsetdf_prefix = ds + '_gp_subset'
        gp_stats_samplesubset.append(get_estimation_stats(gp_subset, 'gp', 'unif_posterior', 'gp_coverage', gp_subsetdf_prefix))

        vbpi_subset = vbpi_df.loc[vbpi_df['exp_samples'] >= sample_min]
        vbpi_subsetdf_prefix = ds + '_vbpi_subset'
        vbpi_stats_samplesubset.append(get_estimation_stats(vbpi_subset, 'vbpi', 'exp_posterior', 'vbpi_coverage', vbpi_subsetdf_prefix))

    gp_full_stats_df = pd.DataFrame(gp_stats_full)
    vbpi_full_stats_df = pd.DataFrame(vbpi_stats_full)
    full_stats_df = gp_full_stats_df.merge(vbpi_full_stats_df, on = 'dataset')
    full_stats_df.to_csv(datapath + "/fulldata_summary_stats.csv", index = False, float_format = '%.5f')

    gp_subset_stats_df = pd.DataFrame(gp_stats_samplesubset)
    vbpi_subset_stats_df = pd.DataFrame(vbpi_stats_samplesubset)
    subset_stats_df = gp_subset_stats_df.merge(vbpi_subset_stats_df, on = 'dataset')
    subset_stats_df.to_csv(datapath + "/subsetdata_summary_stats.csv", index = False, float_format = '%.8f')


    # greater_than = " ".join([">", str(sample_min), "samples"])
    # less_than = " ".join(["<=", str(sample_min), "samples"])

    # full['gtsamplemin'] = np.where(full['unif_samples'] >= sample_min, greater_than, less_than)

    # coverage = full.groupby(['dataset', 'gtsamplemin', 'gp_coverage']).size().agg(
    #     {'count': lambda x: int(x), 'prop': lambda x: x / x.groupby(['dataset', 'gtsamplemin']).sum() * 100}
    #     ).unstack(level=0).reset_index()

    # barplot = (
    #     p9.ggplot(coverage, p9.aes(x = 'dataset', y = 'prop', fill = 'gp_coverage', label = 'prop')) 
    #     + p9.geom_bar(stat=p9.stat_identity)
    #     + p9.geom_text(size = 5, position = p9.position_stack(vjust = 0.5),
    #                 format_string='{:.1f}%')
    #     + p9.facet_wrap('gtsamplemin')
    #     + p9.ylab('PCSP proportion')
    #     + p9.theme(axis_text_x = p9.element_text(size = 6))
    # )
    # barplot.save(datapath + "/coverage.pdf")


def output_stats(datapath, sample_min):
    compile_estimation_stats(datapath, sample_min)
    timing_full = compile_timing_stats(datapath, uniq = False)
    timing_uniq = compile_timing_stats(datapath, uniq = True)

    timing_stats = timing_full.merge(timing_uniq, on = 'dataset')
    timing_stats.to_csv(datapath + '/timing_stats.csv', index = False, float_format = '%.4f')

