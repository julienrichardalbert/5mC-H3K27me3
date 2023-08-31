#!/bin/bash
#SBATCH --ntasks=8                      # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes

# load functions coded in separate files
source ./sra2bw_functions.sh
SCRATCH="/scratch2"

# ---------------- CONFIGURE VARIABLES ----------------- #
FLAG=1540
MIN_MAPQ=10
THREADS=8
THREAD_MEM=4G
BLACKLIST="/data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed"

# ------------------ CALL FUNCTIONS -------------------- #
# takes as input the full path to the fastq file
setupVariables $1
trimReads 5
run_fastQC
run_fastqScreen
alignBowtie2 /data/reference_genomes/mm10/mm10
groomSam
trueStats
bamToBigWigDeeptoolsCPMsmoothKeepDup 1 0
# callPeaks broad
# calculateEnrichment
rename_cleanup
