"""Functions to handle bitset representations of subsplits and GPCSPs."""


def is_rootsplit(gpcsp):
    """Is this "pretty" GPCSP a rootsplit?"""
    return "|" not in gpcsp


def split_pcsp_into_list(pcsp):
    """Split a "pretty" string representation of a PCSP into three components."""
    pcsp_list = pcsp.split("|")
    assert len(pcsp_list) == 3
    return pcsp_list


def parent_split_size(gpcsp):
    """Count the number of taxa in the parent split of a PCSP, or return 1+taxon_count if
    we have a rootsplit."""
    if is_rootsplit(gpcsp):
        return 1 + len(gpcsp)
    # else
    pcsp_list = split_pcsp_into_list(gpcsp)
    return (pcsp_list[0] + pcsp_list[1]).count("1")


def child_subsplit_taxon_set_sizes(gpcsp):
    """Given a "pretty" string representation of a PCSP, find the number of taxa in the
    smallest component of the child subsplit."""
    if is_rootsplit(gpcsp):
        child_component_size_1 = gpcsp.count("0")
        child_component_size_2 = gpcsp.count("1")
    else:
        pcsp_list = split_pcsp_into_list(gpcsp)
        child_component_size_2 = pcsp_list[2].count("1")
        child_component_size_1 = pcsp_list[1].count("1") - child_component_size_2
    return (
        min(child_component_size_1, child_component_size_2),
        max(child_component_size_1, child_component_size_2),
    )


def taxon_set_of_str_bitset(bitset):
    """Make the 1-indexed taxon set from a string bitset."""
    return {index + 1 for index, value in enumerate(bitset) if value == "1"}


def parent_bitset_of_gpcsp(gpcsp_str, sort=False):
    """Given a PCSP string representation, extract the parent subsplit with no bar.
    Return empty string if there is no bar in the input (a rootsplit)."""
    if is_rootsplit(gpcsp_str):
        return ""
    # else
    parent0, parent1, _ = split_pcsp_into_list(gpcsp_str)
    subsplits = [parent0, parent1]
    if sort:
        subsplits.sort()
    return "".join(subsplits)


def rotate_subsplit(bitset_string):
    """Given a bitset abcdef, return defabc."""
    assert len(bitset_string) % 2 == 0
    half_len = len(bitset_string) // 2
    return bitset_string[half_len:] + bitset_string[:half_len]


def pretty_gpcsp_of_string(bitset_string):
    """Insert vertical bars to make a PCSP bitset string easier to read."""
    length = len(bitset_string)
    assert length % 3 == 0
    chunk_length = length // 3
    return "|".join(
        [bitset_string[i : i + chunk_length] for i in range(0, length, chunk_length)]
    )


def indexed_pcsp_of_bitset_pcsp(pcsp_str):
    """Make a string representation of a PCSP described in 1-indexed form."""
    [sister, focal, child1] = [
        taxon_set_of_str_bitset(bitset) for bitset in pcsp_str.split("|")
    ]
    return str(((sister, focal), (set(sorted(focal - child1)), child1)))


def indexed_gpcsp_of_bitset_gpcsp(gpcsp_str):
    """Make a string representation of a pretty GPCSP described in 1-indexed form."""
    if is_rootsplit(gpcsp_str):
        return str(taxon_set_of_str_bitset(gpcsp_str))
    # else
    return indexed_pcsp_of_bitset_pcsp(gpcsp_str)
