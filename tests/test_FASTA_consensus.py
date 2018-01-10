from collections import namedtuple

from scripts import FASTA_consensus


class TestMergeFASTAs:
    def test_with_headers(self):
        sequences = [">seq1", "ACTG",
                     ">seq2", "AATG",
                     ">seq3", "ACAG",
                     ]
        output = {
            0: {'A': 3},
            1: {'C': 2, 'A': 1},
            2: {'T': 2, 'A': 1},
            3: {'G': 3},
            }

        assert FASTA_consensus.merge_FASTAs(sequences) == output

    def test_again(self):
        sequences = ["AATG",
                     "AAAG",
                     "ACAG",
                     ]
        output = {
            0: {'A': 3},
            1: {'C': 1, 'A': 2},
            2: {'T': 1, 'A': 2},
            3: {'G': 3},
            }

        assert FASTA_consensus.merge_FASTAs(sequences) == output

    def test_unequal_input(self):
        sequences = ["AATGA",
                     "AAAG",
                     "ACA",
                     ]
        output = {
            0: {'A': 3},
            1: {'C': 1, 'A': 2},
            2: {'T': 1, 'A': 2},
            3: {'G': 2},
            4: {'A': 1},
            }

        assert FASTA_consensus.merge_FASTAs(sequences) == output


class TestMakeConsensus:
    def test_simple_case(self):
        input_dict = {
            0: {'A': 3},
            1: {'C': 2, 'A': 1},
            2: {'T': 2, 'A': 1},
            3: {'G': 3},
            }
        output = "ACTG"

        assert FASTA_consensus.make_consensus(input_dict) == output

    def test_again(self):

        input_dict = {
            0: {'A': 3},
            1: {'C': 1, 'A': 2},
            2: {'T': 1, 'A': 2},
            3: {'G': 3},
            }
        output = "AAAG"

        assert FASTA_consensus.make_consensus(input_dict) == output

    def test_unequal_input(self):
        input_dict = {
            0: {'A': 3},
            1: {'C': 1, 'A': 2},
            2: {'T': 1, 'A': 2},
            3: {'G': 2},
            4: {'A': 1},
            }
        output = "AAAGA"

        assert FASTA_consensus.make_consensus(input_dict) == output

    def test_gap_not_returned(self):
        input_dict = {
            0: {'A': 3},
            1: {'C': 1, '-': 2},
            2: {'T': 1, 'A': 2},
            3: {'G': 2},
            4: {'A': 1},
            }
        output = "ACAGA"

        assert FASTA_consensus.make_consensus(input_dict) == output

    def test_gap_returned(self):
        input_dict = {
            0: {'A': 3},
            1: {'-': 3},
            2: {'T': 1, 'A':2},
            3: {'G': 2},
            4: {'A': 1},
            }
        output = "A-AGA"

        assert FASTA_consensus.make_consensus(input_dict) == output


frequency = namedtuple('base', 'position A C G T gap depth')

class TestGetBaseFrequency:
    def test_basic(self):
        input_dict = {'C': 2, 'A': 1}

        output = frequency(position=0,
                           A=33.33,
                           C=66.67,
                           G=0,
                           T=0,
                           gap=0,
                           depth=3)

        assert FASTA_consensus.get_base_frequency(input_dict, 0) == output

    def test_gap_merge(self):
        input_dict = {'A': 2, 'C': 3, 'G': 4, 'T': 5, '-': 6, 'N': 7}

        output = frequency(position=1,
                           A=7.41,
                           C=11.11,
                           G=14.81,
                           T=18.52,
                           gap=48.15,
                           depth=27)

        assert FASTA_consensus.get_base_frequency(input_dict, 1) == output

class TestMakeFrequencyMatrix:
    def test_two_positions(self):
        input_dict = {
            0: {'C': 2, 'A': 1},
            1: {'A': 2, 'C': 3, 'G': 4, 'T': 5, '-': 6, 'N': 7},
        }

        freq1 = frequency(position=0,
                          A=33.33,
                          C=66.67,
                          G=0,
                          T=0,
                          gap=0,
                          depth=3)
        freq2 = frequency(position=1,
                          A=7.41,
                          C=11.11,
                          G=14.81,
                          T=18.52,
                          gap=48.15,
                          depth=27)

        output = [freq1, freq2]

        assert FASTA_consensus.make_frequency_matrix(input_dict) == output