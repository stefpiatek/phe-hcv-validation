# HCV pipeline evaluation

Project at Public Health England for STP Elective. 

De novo HCV pipeline used to determine species, quasiaspecies and drug resistance set up to replace Sanger sequencing. As Sanger is not as sensitive as the NGS-based de novo pipeline, this makes an evaluation of determining the sample composition against a known dataset unavilable. 

## Starting dataset

FASTAs were made to simulate samples with different compositions. These have been turned into fastqs using [ART](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3278762/).


## Overview of analysis

- Create a consensus sequence from the FASTAs themselves, and also create a frequency matrix from them (scripts/FASTA_consensus.py, run by process_samples.py)
- Then align the synthetic FASTQs (used BWA) to the conensus sequence and then create a quasibam for the consensus (process_samples.py)
- Compare the frequency differences between the consensus frequency matrix and the consensus quasibam (scripts/frequency-matrix_quaisbam_comparison.Rmd)
- Compare each pipeline quasibam to the consensus frequency matrix (quasibam-pipeline_comparison.Rmd)

## Requirements

- Python 3.5+ required for scripts, have not tested on anything lower. Requirements for pip in pip_requirements.txt
- R 3.2+
    - Ideally I'd use packrat or an anaconda package, but issues with proxy
    - Tidyverse and Biostrings (Biostrings is a bioconductor package)
    - Rstudio if doing single Rmarkdown reports
- Pandoc for knitting Rmarkdown
- For consensus gap, the following in `pipeline-resources/`
    - `hcv.fasta` and associated files
    - `lastz-distrib/bin/lastz` and associated files
    - `lastz_bestref.pl`
    - `lastz_analyser.WITH_REVCOMP.pl`
    - `genome_maker2b.pl`
    - `cons_mv.pl`
    - `N_remover_from_consensus.pl`
    - `majvarcheck2.pl`
- Module files show versions of each tool used

## Running stages

- If running any processing, do this on the cluster
- Then copy `results` folder to local machine and run the report generation
    - Had issues with pandoc and R package installation on the cluster and 
    didn't have time to smooth it out.

#### Default usage

- **Requires** Data in `results/`:
    - Samples following the same YYMMDD_N prefix (N = sample number):
    - FASTA file with each individual synthetic virus named in `{YYMMDD_N}_quasi.fas` format 
    - FASTQ files from derived from the FASTA named in 
     `{YYMMDD_N}_R1.fq` and `{YYMMDD_N}_R2.fq` formats.
- Load modules:
    ```
    module load anaconda/4.2.1_python3
    module load phe/quasi_bam/2-3
    module load samtools # 0.1.19-44428cd
    module load bwa # 0.7.9a-r786
    module load smalt/0.7.6
    module load blast+/2.2.27
    module load vphaser # 2.0
    ```
- Navigate to root directory of the project.
- Run sample processing on the cluster, e.g. 
`python3 process_samples.py 170908`

#### Human removal

- Same as default but run with *--remove-human* flag 
`python3 process_samples.py 171009 --remove-human`


#### Analysing processed data: Rmarkdown reports

- **Requires** the following in the `data` folder:
    - Output from samples to be processed by default or remove-human settings    
    - Quasibam .txt output from pipelines 
    e.g. `171009_12_vicuna_bwa_quasibam.txt`. This format is required.
    - If running another batch of **171009** samples, change the 
    `has_true_positives <- params$date_prefix == "171009"` in the first code chunk
    to be the new date prefix
    - If vphaser has not been run, the vphaser section will be skipped
- Navigate to root directory of the project.
- Run python script with *--reports* flag, e.g. `python3 process_samples.py 170908 --reports`
    - If you want to run a single pipeline only, adding the *--pipeline* flag and specify the pipeline combination. 
    e.g. `python3 process_samples.py 170908 --reports --pipeline vicuna_bwa`
        - This is part of the filename, see the data requirements.
    - You can also run a single pipeline by opening up the
    `scripts/quasibam-pipeline_comparison.Rmd` file in Rstudio, 
    editing the header parameters and then knitting the document in Rstudio.
- Reports will be generated in the `reports` folder, 
look at the html or add the folders and md to gitlab/github for them to be rendered.

#### Vphaser

- **Requires**
    - Samples to be processed by default or remove-human settings
    - Loading of modules listed in default

#### Consensus gap

- **Requires**: Same as default usage
- Navigate to root directory of the project
- Run sample processing on the cluster 
`python3 process_samples.py 170908 --consensus-gap`
- Copy `results/gap_files` to local machine's project
- Open the `scripts/consensus-gap-filling.Rmd` in Rstudio and edit the 
`all_samples <- paste0("170908_", 1:18)` line to fit the correct samples.
    - Then knit in Rstudio


## Example usage

On cluster:

- `python3 process_samples.py 170908`
- `python3 process_samples.py 170908 --vphaser`
- `python3 process_samples.py 171009 --remove-human`

Locally after copying over `results/` contents:

- `python3 process_samples.py 170908 --reports`
- `python3 process_samples.py 171009 --reports --pipeline iva_bwa`


## Testing

Python unit tests were made for FASTA_conensus.py. From the root of the project run:

- `python3 -m pytest tests/` or `pytest tests/`

R unit tests were made for helper_functions. From the root of the project run:

- `Rscript tests/test_helper_functions.R`