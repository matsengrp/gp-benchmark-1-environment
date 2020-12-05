"""our libsbn interface."""

import libsbn


def fit(newick_path, fasta_path, out_csv_path, tol, max_iter):
    inst = libsbn.gp_instance("mmap.dat")
    inst.read_fasta_file(fasta_path)
    inst.read_newick_file(newick_path)
    inst.make_engine()
    inst.print_status()
    inst.estimate_branch_lengths(tol, max_iter)
    inst.estimate_sbn_parameters()
    inst.sbn_parameters_to_csv(out_csv_path)
