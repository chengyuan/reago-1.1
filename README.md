reago v1.1
=====
Changes:
1. Refactored the entire workflow, making code easier to follow.
2. Fixed various bugs.


reago v1.0
=====
an assembly tool for 16S ribosomal RNA recovery from metagenomic data

Dependencies:
1. Infernal 1.1.1
2. Readjoiner 1.2

Runbook:

Input: Paired-end metagenomic reads in FASTA format

Output: 16S genes recovered from metagenomic reads


Step 1: Identify 16S reads.

command: python filter_input.py paired_end_1.fasta paired_end_2.fasta output_dir cm_dir cm_to_use num_of_CPU

example: python filter_input.py sample_1.fasta sample_2.fasta filter_out cm ba 10


Step 2: Assemble 16S reads.

command: python reago.py filename.fasta -l READ_LENGTH

Optional parameters:

-o OVERLAP, default 0.8

-e ERROR_CORRECTION_THRESHOLD, default 0.05

-t TIP_SIZE, default 30

-b PATH_FINDING_PARAMETER, default 10

example: python reago.py make_sample/sample/filtered.fasta sample_out -l 101

**Note:**

REAGO assumes the sequence names of a read pair to be

XXXX.1 &
XXXX.2

The only difference of their names is the last character. 


