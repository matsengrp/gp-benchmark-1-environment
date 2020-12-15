"""our libsbn interface."""

import libsbn
import pandas as pd


def gp_fit(newick_path, fasta_path, out_csv_path, tol, max_iter):
    """Fit an SBN via GP."""
    inst = libsbn.gp_instance("mmap.dat")
    inst.read_fasta_file(fasta_path)
    inst.read_newick_file(newick_path)
    inst.make_engine()
    inst.print_status()
    inst.estimate_branch_lengths(tol, max_iter)
    inst.estimate_sbn_parameters()
    inst.sbn_parameters_to_csv(out_csv_path)


def simple_average(newick_path, out_csv_path):
    """Fit an SBN using simple average."""
    inst = libsbn.rooted_instance("")
    inst.read_newick_file(newick_path)
    inst.process_loaded_trees()
    inst.train_simple_average()
    inst.sbn_parameters_to_csv(out_csv_path)


def tree_probability(newick_path, sbn_parameter_csv, out_csv_path):
    """Calculate tree probabilities."""
    inst = libsbn.rooted_instance("")
    inst.read_newick_file(newick_path)
    inst.process_loaded_trees()
    inst.read_sbn_parameters_from_csv(sbn_parameter_csv)
    probabilities = pd.Series(inst.calculate_sbn_probabilities())
    probabilities.to_csv(out_csv_path, header=False)
