from argparse import ArgumentParser
from collections import defaultdict, namedtuple

"""Simple script to generate consensus sequence from multiple FASTAs

- Multiple FASTAs in a single file, require an output of the most common
  sequence
- merges all FASTAs into a dictionary of base frequency per position,
  then takes most common base per position

Tested with python 3.5.2, using pytest for unit testing
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
        # remove trailing whitespace
        line = line.rstrip()
        for position, base in enumerate(line):
            if base not in sequence_dict[position].keys():
                sequence_dict[position][base] = 1
            else:
                sequence_dict[position][base] += 1
    return sequence_dict


def make_consensus(sequence_dict):
    """Create consesus sequence from sequence_dict

    If maximum base is a gap ('-' or N) and other bases are present,
    use the bases in preference to the gap.

    :param sequence_dict - frequency of each base per sequence postion:
    :return consensus - single string of consensus sequence:
    """

    consensus_list = []

    for position in sequence_dict:
        max_base = max(sequence_dict[position],
                       key=sequence_dict[position].get)
        if max_base not in ['-', 'N']:
            # do not add gaps to consensus
            consensus_list.append(max_base)

    return ''.join(consensus_list)


frequency = namedtuple('base', 'Pos A C G T Gap Depth RefN')
default_frequency = frequency(None, 0, 0, 0, 0, 0, 0, None)


def get_base_frequency(position_dict, position):
    """Make dictionary of base frequency for a single position

    :param position_dict - count of each base for a position:
    :param position - 1-indexed nucleotide position:
    :returns base_frequency - dictionary of base frequencies:
    """
    depth = sum(position_dict.values())
    # combine gaps to single format
    gap = 0
    for gap_notation in ['-', 'N']:
        try:
            gap += position_dict.pop(gap_notation)
        except KeyError:
            continue
    position_dict['Gap'] = gap

    # change counts to frequency by 2 decimal places
    for base, base_count in position_dict.items():
        position_dict[base] = round(100 * base_count / depth, 2)

    position_dict['RefN'] = max(position_dict, key=position_dict.get)
    position_dict['Pos'] = position
    position_dict['Depth'] = depth

    # replace default frequency with position_dict values
    base_frequency = default_frequency._replace(**position_dict)

    return base_frequency


def make_frequency_matrix(sequence_dict):
    """Makes list of base frequencies per position

    - Skips any position where the maximum base is a '-'' or 'N'

    :param sequence_dict - count of each base per position:
    :returns frequency_matrix - list of dictionaries:
    """
    frequency_matrix = []
    base_position = 1
    for position_dict in sequence_dict.values():
        if max(position_dict, key=position_dict.get) in ['-', 'N']:
            # skip any positions where a gap is the dominant base
            continue
        else:
            frequency_matrix.append(
                get_base_frequency(position_dict, base_position)
                )
            base_position += 1

    return frequency_matrix


if __name__ == '__main__':
    # set up argument parser
    parser = ArgumentParser(
        description='Convert FASTAs in one file to a consensus FASTA\n'
                    'Output file made in the same directory as the input')
    parser.add_argument('input_file')
    parser.add_argument('-g', '--gap',
                        help="Insert a gap into consensus sequence")
    parser.add_argument('--gap-sample', default="180212_1",
                        help=("Insert a gap into consensus sequence "
                              "using the given sample name"))
    args = parser.parse_args()
    if args.gap:
        consensus_gap = int(args.gap)
    else:
        consensus_gap = None

    in_file = args.input_file
    consensus_out_file = "{prefix}_consensus.fas".format(
        prefix=in_file.replace('.fas', ''))
    matrix_out_file = "{prefix}_frequency_matrix.txt".format(
        prefix=in_file.replace('.fas', ''))

    assert in_file.endswith(".fas"), "Input file must end with '.fas'"

    # parse input data
    with open(in_file, "r") as sequences:
        sequence_dict = merge_FASTAs(sequences)
        consensus = make_consensus(sequence_dict)
        frequency_matrix = make_frequency_matrix(sequence_dict)
    # write consensus FASTA
    with open(consensus_out_file, "w") as output:
        if consensus_gap is None:
                # No gap added so write entire file
                output.write(">consensus\n")
                output.write(consensus)
        else:
            output.write(">{sample}_quasi_consensus.1\n".format(
                sample=args.gap_sample))
            output.write(consensus[:5000])
            output.write("\n")
            output.write(">{sample}_quasi_consensus.2\n".format(
                sample=args.gap_sample))
            output.write(consensus[5000 + consensus_gap:])
        output.write("\n")
    # write frequency matrix
    with open(matrix_out_file, "w") as output:
        # write header
        output.write('\t'.join(frequency._fields))
        output.write('\n')
        # write data
        for base_frequency in frequency_matrix:
            output.write('\t'.join([str(field)
                                    for field in base_frequency]))
            output.write('\n')
