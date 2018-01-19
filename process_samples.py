import subprocess
from argparse import ArgumentParser
from glob import glob
from os import getcwd

"""
Simple script to run through each sample and take from FASTA to quasibam

At PHE need to load the following modules first

module load anaconda/4.2.1_python3
module load phe/quasi_bam/2-3
module load samtools
module load bwa
"""


directory = getcwd()

parser = ArgumentParser(
    description='Run all processing of FASTAs to create consensus, '
                'frequency matrix and quasibams for comparison')

parser.add_argument('date_prefix', help="Date prefix for samples in YYMMDD")
args = parser.parse_args()
prefix = args.date_prefix

files = glob("{directory}/data/{prefix}_*_quasi.fas".format(
    directory=directory, prefix=prefix))

sample_numbers = [file.split("{}_".format(prefix))[1].split("_")[0]
                  for file in files]


# run all data processing
for sample_number in sample_numbers:
    # create consensus and frequency matrix
    path_prefix = "{directory}/data/{prefix}_{sample_number}".format(
        directory=directory,
        prefix=prefix,
        sample_number=sample_number)

    print("-- Running FASTQ_consensus for sample {sample_number}".format(
        sample_number=sample_number))
    subprocess.run([
        "python3",
        "{directory}/scripts/FASTA_consensus.py".format(
            directory=directory),
        path_prefix + "_quasi.fas"],
        check=True)

    print("-- Running bwa index for sample {sample_number}".format(
        sample_number=sample_number))
    subprocess.run(["bwa", "index",
                    path_prefix + "_quasi_consensus.fas"],
                   check=True)

    print("-- BWA mem for sample {sample_number}".format(
        sample_number=sample_number))

    output_filename = path_prefix + "_quasi.sam"
    with open(output_filename, "w") as output_file:
        subprocess.run(["bwa", "mem",
                        path_prefix + "_quasi_consensus.fas",
                        path_prefix + "_quasi.fas_R1.fq",
                        path_prefix + "_quasi.fas_R2.fq"],
                       stdout=output_file,
                       check=True)

    output_filename = path_prefix + "_quasi.bam"
    with open(output_filename, "w") as output_file:
        print("-- Converting sam to bam for sample {sample_number}".format(
            sample_number=sample_number))
        subprocess.run(["samtools", "view", "-Sb",
                        path_prefix + "_quasi.sam"],
                       stdout=output_file,
                       check=True)

    print("-- Sorting bam for sample {sample_number}".format(
        sample_number=sample_number))
    subprocess.run(["samtools", "sort", "-f",
                    path_prefix + "_quasi.bam",
                    path_prefix + "_quasi_sorted.bam"],
                   check=True)

    print("-- Indexing bam for sample {sample_number}".format(
        sample_number=sample_number))
    subprocess.run(["samtools", "index",
                    path_prefix + "_quasi.bam",
                    path_prefix + "_quasi_sorted.bam"],
                   check=True)

    print("-- Running quasi_bam for sample {sample_number}".format(
        sample_number=sample_number))
    subprocess.run(
        ["quasi_bam",
         # quasi_bam gets path prefix by splitting by ".",
         # so full path can't be given (username contains .)
         "data/{prefix}_{sample_number}_quasi_sorted.bam".format(
             prefix=prefix, sample_number=sample_number),
         "data/{prefix}_{sample_number}_quasi_consensus.fas".format(
             prefix=prefix, sample_number=sample_number),
         "-f 0.001"],
        check=True)
