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

## Testing

Python unit tests were made for FASTA_conensus.py

    `python -m pytest tests/`


## Requirements / versions used

- Python 3.5+ required for scripts, have not tested on anything lower. Requirements for pip in pip_requirements.txt
- samtools 0.1.19
- quasi_bam 2-3
- bwa 0.7.9a

