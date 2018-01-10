from argparse import ArgumentParser
from collections import defaultdict
from sys import version

"""Simple script to generate consensus sequence from multiple FASTAs

- Multiple FASTAs in a single file, require an output of the most common
  sequence
- merges all FASTAs into a dictionary of base frequency per position,
  then takes most common base per position
"""

def merge_FASTAs(sequences):
    """Make dictionary of base counts per nt position
    
    :param sequences - open file or list of FASTA sequences:
    :return sequence_dict - dictionary with keys of each position,
                            the values are a dictionary of each base and their
                            frequency:
    """
    sequence_dict = defaultdict(dict)
    for line in sequences:
        # skip fasta header
        if line.startswith(">"):
            continue
        for position, base in enumerate(line):
            if base not in sequence_dict[position].keys():
                sequence_dict[position][base] = 1
            else:
                sequence_dict[position][base] += 1
    return sequence_dict


def make_consensus(sequence_dict):
    """Create consesus sequence from sequence_dict
    
    :param sequence_dict - frequency of each base per sequence postion:
    :return consensus - single string of consensus sequence:
    """

    consensus = ''.join([max(sequence_dict[base], key=sequence_dict[base].get)
                         for base in sequence_dict])

    return consensus


if __name__ == '__main__':
    # set up argument parser
    parser = ArgumentParser(
        description='Convert multiple FASTAs in one file to a consensus FASTA\n'
                    'Output file made in the same directory as the input')
    parser.add_argument('input_file')
    args = parser.parse_args()

    in_file = args.input_file
    out_file = "{prefix}_consensus.fas".format(
        prefix=in_file.replace('.fas', ''))

    assert in_file.endswith(".fas"), "Input file must end with '.fas'"

    # parse input data
    with open(in_file, "r") as sequences:
        sequence_dict = merge_FASTAs(sequences)
        consensus = make_consensus(sequence_dict)
    # write consensus FASTA
    with open(out_file, "w") as output:
        output.write(">consensus\n")
        output.write(consensus)
