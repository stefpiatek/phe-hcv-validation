from collections import namedtuple

from scripts import FASTA_consensus


class TestMergeFASTAs:
    def test_with_headers(self):
        sequences = [">seq1", "ACTG",
                     ">seq2", "AATG",
                     ">seq3", "ACAG",
                     ]
        output = {
            1: {'A': 3},
            2: {'C': 2, 'A': 1},
            3: {'T': 2, 'A': 1},
            4: {'G': 3},
            }

        assert FASTA_consensus.merge_FASTAs(sequences) == output

    def test_again(self):
        sequences = ["AATG",
                     "AAAG",
                     "ACAG",
                     ]
        output = {
            1: {'A': 3},
            2: {'C': 1, 'A': 2},
            3: {'T': 1, 'A': 2},
            4: {'G': 3},
            }

        assert FASTA_consensus.merge_FASTAs(sequences) == output

    def test_unequal_input(self):
        sequences = ["AATGA",
                     "AAAG",
                     "ACA",
                     ]
        output = {
            1: {'A': 3},
            2: {'C': 1, 'A': 2},
            3: {'T': 1, 'A': 2},
            4: {'G': 2},
            5: {'A': 1},
            }

        assert FASTA_consensus.merge_FASTAs(sequences) == output


class TestMakeConsensus:
    def test_simple_case(self):
        input_dict = {
            1: {'A': 3},
            2: {'C': 2, 'A': 1},
            3: {'T': 2, 'A': 1},
            4: {'G': 3},
            }
        output = "ACTG"

        assert FASTA_consensus.make_consensus(input_dict) == output

    def test_again(self):

        input_dict = {
            1: {'A': 3},
            2: {'C': 1, 'A': 2},
            3: {'T': 1, 'A': 2},
            4: {'G': 3},
            }
        output = "AAAG"

        assert FASTA_consensus.make_consensus(input_dict) == output

    def test_unequal_input(self):
        input_dict = {
            1: {'A': 3},
            2: {'C': 1, 'A': 2},
            3: {'T': 1, 'A': 2},
            4: {'G': 2},
            5: {'A': 1},
            }
        output = "AAAGA"

        assert FASTA_consensus.make_consensus(input_dict) == output

    def test_gap_not_returned(self):
        input_dict = {
            1: {'A': 3},
            2: {'C': 1, '-': 2},
            3: {'T': 1, 'A': 2},
            4: {'G': 2},
            5: {'A': 1},
            }
        output = "ACAGA"

        assert FASTA_consensus.make_consensus(input_dict) == output

    def test_gap_returned(self):
        input_dict = {
            1: {'A': 3},
            2: {'-': 3},
            3: {'T': 1, 'A': 2},
            4: {'G': 2},
            5: {'A': 1},
            }
        output = "A-AGA"

        assert FASTA_consensus.make_consensus(input_dict) == output


frequency = namedtuple('base', 'Pos A C G T Gap Depth')


class TestGetBaseFrequency:
    def test_basic(self):
        input_dict = {'C': 2, 'A': 1}

        output = frequency(Pos=0,
                           A=33.33,
                           C=66.67,
                           G=0,
                           T=0,
                           Gap=0,
                           Depth=3)

        assert FASTA_consensus.get_base_frequency(input_dict, 0) == output

    def test_gap_merge(self):
        input_dict = {'A': 2, 'C': 3, 'G': 4, 'T': 5, '-': 6, 'N': 7}

        output = frequency(Pos=1,
                           A=7.41,
                           C=11.11,
                           G=14.81,
                           T=18.52,
                           Gap=48.15,
                           Depth=27)

        assert FASTA_consensus.get_base_frequency(input_dict, 1) == output


class TestMakeFrequencyMatrix:
    def test_two_positions(self):
        input_dict = {
            0: {'C': 2, 'A': 1},
            1: {'A': 2, 'C': 3, 'G': 4, 'T': 5, '-': 6, 'N': 7},
        }

        freq1 = frequency(Pos=0,
                          A=33.33,
                          C=66.67,
                          G=0,
                          T=0,
                          Gap=0,
                          Depth=3)
        freq2 = frequency(Pos=1,
                          A=7.41,
                          C=11.11,
                          G=14.81,
                          T=18.52,
                          Gap=48.15,
                          Depth=27)

        output = [freq1, freq2]

        assert FASTA_consensus.make_frequency_matrix(input_dict) == output
