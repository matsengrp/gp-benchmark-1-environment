"""our libsbn interface."""

import glob
import libsbn
import numpy as np
import os
import pandas as pd
from scipy.special import logsumexp, softmax


def gp_fit(newick_path, fasta_path, out_csv_prefix, tol, max_iter):
    """Fit an SBN via GP."""
    inst = libsbn.gp_instance("mmap.dat")
    inst.read_fasta_file(fasta_path)
    inst.read_newick_file(newick_path)
    inst.make_engine()
    inst.print_status()
    inst.estimate_branch_lengths(tol, max_iter)
    inst.estimate_sbn_parameters()
    inst.sbn_parameters_to_csv(out_csv_prefix + ".sbn.csv")
    inst.branch_lengths_to_csv(out_csv_prefix + ".bl.csv")


def simple_average(newick_path, out_csv_prefix):
    """Fit an SBN using simple average."""
    inst = libsbn.rooted_instance("")
    inst.read_newick_file(newick_path)
    inst.process_loaded_trees()
    inst.train_simple_average()
    inst.sbn_parameters_to_csv(out_csv_prefix + ".csv")
    inst.unconditional_subsplit_probabilities_to_csv(out_csv_prefix + ".subsplit.csv")


def tree_probability(newick_path, sbn_parameter_csv, out_csv_path):
    """Calculate tree probabilities."""
    inst = libsbn.rooted_instance("")
    inst.read_newick_file(newick_path)
    inst.process_loaded_trees()
    inst.read_sbn_parameters_from_csv(sbn_parameter_csv)
    probabilities = pd.Series(inst.calculate_sbn_probabilities())
    probabilities.to_csv(out_csv_path, header=False, index=False)


def marginal_across_newick_trees(newick_path, fasta_path):
    """Directly (i.e. using an average) estimate the marginal log likelihood for trees
    supplied in a file."""
    simple_specification = libsbn.PhyloModelSpecification(
        substitution="JC69", site="constant", clock="none"
    )
    inst = libsbn.unrooted_instance("")
    inst.read_fasta_file(fasta_path)
    inst.read_newick_file(newick_path)
    inst.prepare_for_phylo_likelihood(simple_specification, 2, [])
    likelihoods = np.array(inst.log_likelihoods())
    return logsumexp(likelihoods) - np.log(len(likelihoods))


def tree_marginal(newick_glob, fasta_path, out_csv_path):
    """Directly estimate the marginal log likelihood for trees supplied in a file for
    each file in the supplied Newick path glob."""
    paths = glob.glob(newick_glob)
    marginals = [
        marginal_across_newick_trees(newick_path, fasta_path) for newick_path in paths
    ]
    prefixes = [os.path.basename(path).split(".")[0] for path in paths]

    df = pd.DataFrame({"prefixes": prefixes, "marginals": marginals})
    df.sort_values(df.columns.values.tolist(), inplace=True)
    df.to_csv(out_csv_path, index=False)


def export_trees_with_subsplits(newick_path, fasta_path, pcsp_csv_path, tol, max_iter):
    """Fit a GP with the given sequences and trees, then for every given PCSP export the
    trees from `newick_path` with the GP branch lengths to `_ignore/trees/...`.

    We supply the PCSPs as the first column of the given CSV of PCSPs (typically the SBN
    parameters from a run of GP).
    """
    inst = libsbn.gp_instance("mmap.dat")
    inst.read_fasta_file(fasta_path)
    inst.read_newick_file(newick_path)
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
