from scripts import FASTA_consensus



class TestMergeFASTAs:
    def test_with_headers(self):
        sequences = [">seq1", "ACTG",
                     ">seq2", "AATG",
                     ">seq3", "ACAG",
                     ]
        output = {
            0: {'A': 3},
            1: {'C': 2, 'A':1},
            2: {'T': 2, 'A':1},
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
            1: {'C': 1, 'A':2},
            2: {'T': 1, 'A':2},
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
            1: {'C': 1, 'A':2},
            2: {'T': 1, 'A':2},
            3: {'G': 2},
            4: {'A': 1},
            }

        assert FASTA_consensus.merge_FASTAs(sequences) == output


class TestMakeConsensus:
    def test_simple_case(self):
        input_dict = {
            0: {'A': 3},
            1: {'C': 2, 'A':1},
            2: {'T': 2, 'A':1},
            3: {'G': 3},
            }
        output = "ACTG"

        assert FASTA_consensus.make_consensus(input_dict) == output

    def test_again(self):

        input_dict = {
            0: {'A': 3},
            1: {'C': 1, 'A':2},
            2: {'T': 1, 'A':2},
            3: {'G': 3},
            }
        output = "AAAG"

        assert FASTA_consensus.make_consensus(input_dict) == output

    def test_unequal_input(self):
        input_dict = {
            0: {'A': 3},
            1: {'C': 1, 'A':2},
            2: {'T': 1, 'A':2},
            3: {'G': 2},
            4: {'A': 1},
            }
        output = "AAAGA"

        assert FASTA_consensus.make_consensus(input_dict) == output
