import subprocess
from argparse import ArgumentParser
from glob import glob
from os import getcwd, makedirs
from os.path import dirname, exists
from shutil import rmtree

"""
Simple script to run through each sample and take from FASTA to quasibam

At PHE need to load the following modules first

module load anaconda/4.2.1_python3
module load phe/quasi_bam/2-3
module load samtools
module load bwa
module load smalt/0.7.6
module load blast+/2.2.27
module load vphaser
"""

directory = getcwd()

parser = ArgumentParser(
    description='Run all processing of FASTAs to create consensus, '
                'frequency matrix and quasibams for comparison')
parser.add_argument('date_prefix', help="Date prefix for samples in YYMMDD")
parser.add_argument('--remove-human',
                    help='Carry out pipeline step of removing human reads',
                    action='store_true')
parser.add_argument('-r', '--reports',
                    help='Render Rmd reports only',
                    action='store_true')
parser.add_argument('-fm', '--fastq-middle', default="_",
                    help=("string between 'YYMMDD_sample_number' and '.fq'. "
                          " e.g. '_quasi.fas_'"))
parser.add_argument('--consensus-gap', action='store_true',
                    help=("Carry out steps to test whether gaps in consensus "
                          "contigs are merged from pipeline steps"))
parser.add_argument('--vphaser', action='store_true',
                    help=("Run vphaser the sample set"))

args = parser.parse_args()
prefix = args.date_prefix
fastq_middle = args.fastq_middle

# Rmd reports only --

if args.reports:
    cmd = ("Rscript -e \"rmarkdown::render("
           "'scripts/frequency-matrix_comparison.Rmd', "
           "params = list(date_prefix = '{prefix}'), "
           "'html_document', "
           "'../reports/{prefix}_frequency-matrix_comparison.html')\""
           )
    subprocess.run(cmd.format(prefix=prefix),
                   shell=True, check=True)
    for pipeline in ["vicuna_bwa", "vicuna_smalt", "iva_bwa", "iva_smalt"]:

        cmd = ("Rscript -e \"rmarkdown::render("
               "'scripts/quasibam-pipeline_comparison.Rmd', "
               "params = list(pipeline = '{pipeline}', "
               "date_prefix = '{prefix}'), "
               "'html_document', "
               "'../reports/{prefix}_{output_pipeline}_report.html')\""
               )
        subprocess.run(
            cmd.format(prefix=prefix,
                       pipeline=pipeline,
                       output_pipeline=pipeline.replace("_", "-")),
            shell=True, check=True)
    exit()

# Processing of files from fasta + fq onwards --

files = glob("{directory}/data/{prefix}_*_quasi.fas".format(
    directory=directory, prefix=prefix))

sample_numbers = [file.split("{}_".format(prefix))[1].split("_")[0]
                  for file in files]

if args.consensus_gap:
    # Pretty messy conversion of bash script so tucked away in
    # another file
    for sample_number in range(1, 19):
        subprocess.run([
            "python3",
            "{directory}/scripts/consensus_gap.py".format(
                directory=directory),
            "{prefix}_{sample_number}".format(
                prefix=prefix,
                sample_number=sample_number)],
            check=True)
    exit()

if args.vphaser:
    for sample_number in sample_numbers:
        print("-- Running vphaser2 for sample {}".format(sample_number))
        vphaser_dir = ("{directory}/data/vphaser/{prefix}_{sample_number}"
                       ).format(
                directory=directory,
                prefix=prefix,
                sample_number=sample_number)
        if not exists(vphaser_dir):
            makedirs(vphaser_dir)
        try:
            subprocess.run([
                "vphaser",
                "-i",
                ("{directory}/data/{prefix}_{sample_number}_quasi_sorted.bam"
                 ).format(directory=directory,
                          prefix=prefix,
                          sample_number=sample_number),
                "-o", vphaser_dir],
                check=True)
        except subprocess.CalledProcessError:
            print("Sample {} failed...\nMoving on...".format(sample_number))
            rmtree(vphaser_dir)
    exit()


# run all data processing
for sample_number in sample_numbers:
    # create consensus and frequency matrix
    sample_prefix = "{directory}/data/{prefix}_{sample_number}".format(
        directory=directory,
        prefix=prefix,
        sample_number=sample_number)

    print("-- Running FASTQ_consensus for sample {sample_number}".format(
        sample_number=sample_number))

    subprocess.run([
        "python3",
        "{directory}/scripts/FASTA_consensus.py".format(
            directory=directory),
        sample_prefix + "_quasi.fas"],
        check=True)

    print("-- Running bwa index for sample {sample_number}".format(
        sample_number=sample_number))
    subprocess.run(["bwa", "index",
                    sample_prefix + "_quasi_consensus.fas"],
                   check=True)

    if args.remove_human:
        # Copied human filtering out from basic pipeline script
        # Kept the same to avoid changing the outcome of filtering
        print("-- Filtering human reads for sample {sample_number}".format(
            sample_number=sample_number))
        # in awk command,  double { to avoid string formatting
        cmd = ("smalt map -x -y 0.5 -i 500 -n 8 "
               "{directory}/pipeline-resources/hg38/hg38_hcv_k15_s3 "
               "{sample_prefix}_R1.fq {sample_prefix}_R2.fq "
               "| awk '{{if ($3 !~ /^chr/ && $7 !~ /^chr/) print $0}}' "
               "> {sample_prefix}_filtered.sam")
        subprocess.run(
            cmd.format(
                directory=directory,
                sample_prefix=sample_prefix),
            shell=True, check=True)

        print("-- convert filtered sam to fastqs for # {sample_number}".format(
            sample_number=sample_number))
        cmd = ("samtools view -bShf 64 {sample_prefix}_filtered.sam"
               "| samtools bam2fq - > "
               "{sample_prefix}{fastq_middle}R1_filtered.fq")
        subprocess.run(
            cmd.format(sample_prefix=sample_prefix, fastq_middle=fastq_middle),
            shell=True, check=True)
        cmd = ("samtools view -bShf 128 {sample_prefix}_filtered.sam"
               "| samtools bam2fq - > "
               "{sample_prefix}{fastq_middle}R2_filtered.fq")
        subprocess.run(
            cmd.format(sample_prefix=sample_prefix, fastq_middle=fastq_middle),
            shell=True, check=True)
        fastq_suffix = "_filtered.fq"
    else:
        fastq_suffix = ".fq"

    print("-- BWA mem for sample {sample_number}".format(
        sample_number=sample_number))

    output_filename = sample_prefix + "_quasi.sam"
    with open(output_filename, "w") as output_file:
        subprocess.run(["bwa", "mem",
                        sample_prefix + "_quasi_consensus.fas",
                        sample_prefix + fastq_middle + "R1" + fastq_suffix,
                        sample_prefix + fastq_middle + "R2" + fastq_suffix],
                       stdout=output_file,
                       check=True)

    output_filename = sample_prefix + "_quasi.bam"
    with open(output_filename, "w") as output_file:
        print("-- Converting sam to bam for sample {sample_number}".format(
            sample_number=sample_number))
        subprocess.run(["samtools", "view", "-Sb",
                        sample_prefix + "_quasi.sam"],
                       stdout=output_file,
                       check=True)

    print("-- Sorting bam for sample {sample_number}".format(
        sample_number=sample_number))
    subprocess.run(["samtools", "sort", "-f",
                    sample_prefix + "_quasi.bam",
                    sample_prefix + "_quasi_sorted.bam"],
                   check=True)

    print("-- Indexing bam for sample {sample_number}".format(
        sample_number=sample_number))
    subprocess.run(["samtools", "index",
                    sample_prefix + "_quasi_sorted.bam"],
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
