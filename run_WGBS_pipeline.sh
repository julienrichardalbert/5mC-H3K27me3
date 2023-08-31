#!/bin/bash

# usage
# run_DNAme_mm10.sh <full path to file>.NO.EXTENSION
# example:
# extension MUST BE "_R1.fastq.gz"

# NOTE: USES 80% OF AVAILABLE RAM FOR METHYLATION EXTRACTOR STEP

# load functions coded in separate files
source ./sra2bw_functions.sh

# ---------------- CONFIGURE VARIABLES ----------------- #
FLAG=1540
MIN_MAPQ=1
THREADS=12
THREAD_MEM=3G
BLACKLIST="/data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed"
# ------------------ CALL FUNCTIONS -------------------- #

# takes as input the full path to the fastq file
setupVariables $1
trimReads 1
run_fastQC
# run_fastqScreen_convertedDNA

alignBismark /data/reference_genomes/mm10/
collectBismarkStats
mergeTwoStrandMethylation
convertMethylationToBigWig 1 /data/reference_genomes/mm10/mm10.sizes
convertMethylationToBigWig 5 /data/reference_genomes/mm10/mm10.sizes

cleanupBismark
