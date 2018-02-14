import subprocess
from os import getcwd, makedirs
from os.path import dirname, exists


def align_and_pileup(align_suffix, bam_suffix):
    """Carry out alignment steps to mpileup
    This is done three times throughout script so avoids some repeats
    :param align_suffix: str, e.g. "_quasi_consensus.fas"
    :param bam_suffix: str, e.g. "_contigs.bam"

    Note: may need to change this to BWA to test this condition?
    """
    subprocess.run(
        ["smalt", "index", "-k", "15", "-s", "3",
         sample_prefix + bam_suffix.replace(".bam", ".k15_s3"),
         sample_prefix + align_suffix],
        check=True)

    subprocess.run(
        ["smalt", "map", "-x", "-y", "0.5", "-i", "500",
         "-n", "8", "-f", "sam",
         "-o", sample_prefix + bam_suffix.replace(".bam", ".sam"),
         sample_prefix + bam_suffix.replace(".bam", ".k15_s3"),
         sample_prefix + "_quasi_R1.fq",
         sample_prefix + "_quasi_R2.fq"],
        check=True)

    output_filename = sample_prefix + bam_suffix
    with open(output_filename, "w") as output_file:
        subprocess.run(["samtools", "view", "-Sb",
                        sample_prefix + bam_suffix.replace(".bam", ".sam")],
                       stdout=output_file,
                       check=True)

    subprocess.run(["samtools", "sort", "-f",
                    sample_prefix + bam_suffix,
                    sample_prefix + bam_suffix.replace(".bam", "sorted.bam")],
                   check=True)

    subprocess.run(["samtools", "index",
                    sample_prefix + bam_suffix.replace(".bam", "sorted.bam")],
                   check=True)

    output_filename = sample_prefix + bam_suffix.replace(".bam", ".mpileup")
    with open(output_filename, "w") as output_file:
        subprocess.run(
            ["samtools", "mpileup",
             "-f", sample_prefix + align_suffix,
             "-d", "1000000",
             sample_prefix + bam_suffix.replace(".bam", "sorted.bam")],
            check=True,
            stdout=output_file)

    cmd = ("awk {{'print $1\"\t\"$2\"\t\"$3\"\t\"$4'}} "
           "{sample_prefix}{pile_suff} > "
           "{sample_prefix}{fake_pile}")

    subprocess.run(
        cmd.format(
            sample_prefix=sample_prefix,
            pile_suff=bam_suffix.replace(".bam", ".mpileup"),
            fake_pile=bam_suffix.replace(".bam", ".fake_mpileup")),
        shell=True, check=True)


directory = getcwd()
prefix = "180212"

gap_folder = dirname("{directory}/data/gap_files/".format(
        directory=directory))
if not exists(gap_folder):
    makedirs(gap_folder)

sample_in = "{directory}/data/170908_1_quasi.fas".format(
    directory=directory)


# iterate through
for sample_number, gap in zip([1, 2, 3, 4, 5, 6, 7],
                              [0, 50, 100, 200, 400, 800, 1600]):
    sample_folder = dirname(
        "{directory}/data/gap_files/{prefix}_{sample_number}/".format(
            directory=directory,
            prefix=prefix,
            sample_number=sample_number))
    if not exists(sample_folder):
        makedirs(sample_folder)

    sample_prefix = (
        "{directory}/data/gap_files/{prefix}_{sample_number}/"
        "{prefix}_{sample_number}"
        ).format(directory=directory,
                 prefix=prefix,
                 sample_number=sample_number)
    resource_prefix = ("{directory}/pipeline-resources/"
                       ).format(directory=directory)

    subprocess.run(["cp", sample_in, sample_prefix + "_quasi.fas"],
                   check=True)

    subprocess.run(["cp", sample_in + "_R1.fq",
                    sample_prefix + "_quasi_R1.fq"],
                   check=True)
    subprocess.run(["cp", sample_in + "_R2.fq",
                    sample_prefix + "_quasi_R2.fq"],
                   check=True)

    print("-- Running FASTQ_consensus for sample {sample_number}".format(
        sample_number=sample_number))

    subprocess.run([
        "python3",
        "{directory}/scripts/FASTA_consensus.py".format(
            directory=directory),
        sample_prefix + "_quasi.fas",
        "-g {gap}".format(gap=gap),
        "--gap-sample", "180212_{sample_number}".format(
            sample_number=sample_number)],
        check=True)

    print("-- Running lastz for sample {sample_number}".format(
        sample_number=sample_number))

    output_filename = sample_prefix + "_contigs.lastz"
    with open(output_filename, "w") as output_file:
        subprocess.run(
            [resource_prefix + "lastz-distrib/bin/lastz",
             sample_prefix + "_quasi_consensus.fas[multiple]",
             resource_prefix + "hcv.fasta",
             "--ambiguous=iupac",
             "--format=GENERAL"],
            stdout=output_file,
            check=True)

    print("-- Analyzing lastz for sample {sample_number}".format(
            sample_number=sample_number))
    subprocess.run(
        ["perl", "-s",
         resource_prefix + "lastz_bestref.pl",
         "-contig_lastz=" + sample_prefix + "_contigs.lastz",
         "-blastdb=" + resource_prefix + "hcv.fasta",
         "-best_ref_fasta=" + sample_prefix + "_ref.fas",
         "-lastz_best_hit_log=" + sample_prefix + "_best_ref.log"],
        check=True)

    print("-- Comparing contigs and best ref: {sample_number}".format(
            sample_number=sample_number))
    output_filename = sample_prefix + "_contig-vs-bestref.lav"
    with open(output_filename, "w") as output_file:
        subprocess.run(
            [resource_prefix + "lastz-distrib/bin/lastz",
             sample_prefix + "_ref.fas",
             sample_prefix + "_quasi_consensus.fas",
             "--ambiguous=iupac"],
            stdout=output_file,
            check=True)

    print("-- Final lastz analysis for sample {sample_number}".format(
            sample_number=sample_number))
    subprocess.run(
        ["perl", "-w", "-s",
         resource_prefix + "lastz_analyser.WITH_REVCOMP.pl",
         "-reference_fasta_file=" + sample_prefix + "_ref.fas",
         "-sample_fasta_file=" + sample_prefix + "_quasi_consensus.fas",
         "-lastz_results_file=" + sample_prefix + "_contig-vs-bestref.lav",
         "-cutoff=50000", "-with_revcomp=yes",
         "-output=" + sample_prefix + "_lastz_analysed_file",
         "-log_file=" + sample_prefix + "_lastz_analyser.log"],
        check=True)

    align_and_pileup(align_suffix="_quasi_consensus.fas",
                     bam_suffix="_contigs.bam")

    print("-- Running genome maker for {sample_number}".format(
        sample_number=sample_number))
    subprocess.run(
        ["perl", "-w", "-s",
         resource_prefix + "genome_maker2b.pl",
         "-sample_pileup_file=" + sample_prefix + "_contigs.mpileup",
         "-contigs=" + sample_prefix + "_quasi_consensus.fas",
         "-reference_mapped_consensus=" + sample_prefix + "_ref.fas",
         "-lastz_analysed_file=" + sample_prefix + "_lastz_analysed_file",
         "-ref_correct_start=0", "-ref_correct_stop=20000",
         "-output=" + sample_prefix + "_genome.fas",
         "-log_file=" + sample_prefix + "_genome_maker_log"],
        check=True)

    print(("-- Iteration 1: consensus from draft genome: {sample_number}"
           ).format(sample_number=sample_number))

    align_and_pileup(align_suffix="_genome.fas",
                     bam_suffix="_genome.bam")

    print("-- Iteration 1: running Cons_mv for {sample_number}".format(
        sample_number=sample_number))
    subprocess.run(
        ["perl", "-w", "-s",
         resource_prefix + "cons_mv.pl",
         "-reference_fasta=" + sample_prefix + "_genome.fas",
         "-mpileup=" + sample_prefix + "_genome.mpileup",
         "-mv_freq_cutoff=0.01",
         "-mv_overall_depth_cutoff=100",
         "-mv_variant_depth_cutoff=20",
         "-cons_depth_cutoff=80",
         "-sliding_window_size=300",
         "-consensus_out=" + sample_prefix + "_consensus1_preNcut.fas",
         "-mv_out=" + sample_prefix + "_genome.fas.mv",
         "-base_freq_out=" + sample_prefix + "_genome.fas.basefreqs.tsv"],
        check=True)

    output_filename = sample_prefix + "_consensus1.fas"
    with open(output_filename, "w") as output_file:
        subprocess.run(
            ["perl", "-w", "-s",
             resource_prefix + "N_remover_from_consensus.pl",
             "-cutoff=46",
             sample_prefix + "_consensus1_preNcut.fas"],
            check=True,
            stdout=output_file)

    print(("-- Iteration 2: consensus from draft genome: {sample_number}"
           ).format(sample_number=sample_number))

    align_and_pileup(align_suffix="_consensus1.fas",
                     bam_suffix="_consensus1.bam")

    print("-- Iteration 2: running Cons_mv for {sample_number}".format(
        sample_number=sample_number))
    subprocess.run(
        ["perl", "-w", "-s",
         resource_prefix + "cons_mv.pl",
         "-reference_fasta=" + sample_prefix + "_consensus1.fas",
         "-mpileup=" + sample_prefix + "_consensus1.mpileup",
         "-mv_freq_cutoff=0.01",
         "-mv_overall_depth_cutoff=100",
         "-mv_variant_depth_cutoff=20",
         "-cons_depth_cutoff=80",
         "-sliding_window_size=300",
         "-consensus_out=" + sample_prefix + "_consensus2_preNcut.fas",
         "-mv_out=" + sample_prefix + "_consensus1.fas.mv",
         "-base_freq_out=" + sample_prefix + "_consensus1.fas.basefreqs.tsv"],
        check=True)

    output_filename = sample_prefix + "_consensus2.fasta"
    with open(output_filename, "w") as output_file:
        subprocess.run(
            ["perl", "-w", "-s",
             resource_prefix + "N_remover_from_consensus.pl",
             "-cutoff=46",
             sample_prefix + "_consensus2_preNcut.fas"],
            check=True,
            stdout=output_file)

    print("-- Running majvarcheck2 for sample {sample_number}".format(
        sample_number=sample_number))
    subprocess.run(
        ["perl", "-w", "-s",
         resource_prefix + "majvarcheck2.pl",
         "-mvpath=" + sample_prefix + "_consensus2.fasta.mv",
         "-basefreq=" + sample_prefix + "_consensus1.fas.basefreqs.tsv",
         "-fwdreads=" + sample_prefix + "_quasi_R1.fq",
         "-revreads=" + sample_prefix + "_quasi_R2.fq"],
        check=True)
