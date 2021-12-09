"""our bito interface."""

import glob
import bito
import numpy as np
import os
import pandas as pd
from scipy.special import logsumexp, softmax


def gp_fit(newick_path, fasta_path, out_csv_prefix, tol, max_iter, use_gradients, mmap_path):
    """Fit an SBN via GP."""
    inst = bito.gp_instance(mmap_path)
    inst.read_fasta_file(fasta_path)
    inst.read_newick_file(newick_path)
    inst.use_gradient_optimization(use_gradients)
    inst.make_engine()
    inst.print_status()
    inst.estimate_branch_lengths(tol, max_iter)
    inst.calculate_hybrid_marginals()
    inst.estimate_sbn_parameters() 
    inst.sbn_parameters_to_csv(out_csv_prefix + ".sbn.csv")
    inst.branch_lengths_to_csv(out_csv_prefix + ".bl.csv")
    inst.sbn_prior_to_csv(out_csv_prefix + ".prior.csv")
    inst.per_gpcsp_llhs_to_csv(out_csv_prefix + ".perpcsp.csv")
    inst.per_gpcsp_bls_from_opt_to_csv(out_csv_prefix + ".perpcsp_bl_from_opt.csv")
    inst.per_gpcsp_llhs_from_opt_to_csv(out_csv_prefix + ".perpcsp_llh_from_opt.csv")


def pcsp_likelihood_surface(newick_path, fasta_path, out_csv_prefix, tol, max_iter, use_gradients, steps, mmap_path):
    """Get the per PCSP log likelihood surfaces and then track perturbations from llh change"""
    inst = bito.gp_instance(mmap_path)
    inst.read_fasta_file(fasta_path)
    inst.read_newick_file(newick_path)
    inst.use_gradient_optimization(use_gradients)
    inst.make_engine()
    inst.print_status()
    inst.hot_start_branch_lengths()

    if steps > 0 & use_gradients == 0:
        inst.scan_pcsp_likelihoods(steps)
        inst.per_gpcsp_llh_surfaces_to_csv(out_csv_prefix + ".perpcsp_llh_surface.csv")
    
    # Now fiddle/track
    inst.hot_start_branch_lengths();
    inst.track_values_from_optimization();
    inst.full_dag_optim_values_to_csv(out_csv_prefix + ".tracked_bl_correction.csv")
    inst.optim_path_bl_to_csv(out_csv_prefix + ".optim_path_bl.csv")
    inst.optim_path_llh_to_csv(out_csv_prefix + ".optim_path_llh.csv")


def simple_average(newick_path, out_csv_prefix):
    """Fit an SBN using simple average."""
    inst = bito.rooted_instance("")
    inst.read_newick_file(newick_path)
    inst.process_loaded_trees()
    inst.train_simple_average()
    inst.sbn_parameters_to_csv(out_csv_prefix + ".csv")
    inst.unconditional_subsplit_probabilities_to_csv(out_csv_prefix + ".subsplit.csv")


def tree_probability(newick_path, sbn_parameter_csv, out_csv_path):
    """Calculate tree probabilities."""
    inst = bito.rooted_instance("")
    inst.read_newick_file(newick_path)
    inst.process_loaded_trees()
    inst.read_sbn_parameters_from_csv(sbn_parameter_csv)
    probabilities = pd.Series(inst.calculate_sbn_probabilities())
    probabilities.to_csv(out_csv_path, header=False, index=False)


def marginal_across_newick_trees(newick_path, fasta_path):
    """Directly (i.e. using an average) estimate the marginal log likelihood for trees
    supplied in a file."""
    simple_specification = bito.PhyloModelSpecification(
        substitution="JC69", site="constant", clock="none"
    )
    inst = bito.unrooted_instance("")
    inst.read_fasta_file(fasta_path)
    inst.read_newick_file(newick_path)
    inst.prepare_for_phylo_likelihood(simple_specification, 2, [])
    likelihoods = np.array(inst.log_likelihoods())
    return logsumexp(likelihoods) - np.log(len(likelihoods))


def tree_marginal(newick_glob, fasta_path, out_csv_path):
    """Directly estimate the marginal log likelihood for trees supplied in a file for
    each file in the supplied Newick path glob."""
    paths = glob.glob(newick_glob)
    assert paths
    marginals = [
        marginal_across_newick_trees(newick_path, fasta_path) for newick_path in paths
    ]
    prefixes = [os.path.basename(path).split(".")[0] for path in paths]

    df = pd.DataFrame({"gpcsp": prefixes, "marginal": marginals})
    df.sort_values(df.columns.values.tolist(), inplace=True)
    df.to_csv(out_csv_path, index=False)


def export_trees_with_subsplits(newick_path, fasta_path, pcsp_csv_path, tol, max_iter, use_gradients):
    """Fit a GP with the given sequences and trees, then for every given PCSP export the
    trees from `newick_path` with the GP branch lengths to `_ignore/trees/...`.

    We supply the PCSPs as the first column of the given CSV of PCSPs (typically the SBN
    parameters from a run of GP).
    """
    inst = bito.gp_instance("mmap.dat")
    inst.read_fasta_file(fasta_path)
    inst.read_newick_file(newick_path)
    inst.use_gradient_optimization(use_gradients)
    inst.make_engine()
    inst.print_status()
    inst.estimate_branch_lengths(tol, max_iter)

    sbn_df = pd.read_csv(pcsp_csv_path)

    for gpcsp in sbn_df.iloc[:, 0]:
        if "|" in gpcsp:
            pcsp_collapsed = gpcsp.replace("|", "")
            inst.export_trees_with_a_pcsp(
                pcsp_collapsed, f"_ignore/trees/{pcsp_collapsed}.nwk"
            )
