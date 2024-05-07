# Code used to make revisions for the DNA methylation - PRC2 paper

# Reviewer 1 comment 1: peak calling sucks and I agree. Let's try improving it anyway
# H3K27me3
macs2 callpeak
    --broad \
    -g mm \
    --keep-dup all \
    --format BAM \
    -t  ESC_serum_TKO_H3K27me3_ChIP_Cooper2014_rep1_ERR476767_trimV5_mm10.bam ESC_serum_TKO_H3K27me3_ChIP_Cooper2014_rep2_ERR476744_trimV5_mm10.bam \
    -c ESC_serum_TKO_input_ChIP_Cooper2014_rep1_ERR476766_ERR476752_trimV5_mm10.bam  ESC_serum_TKO_input_ChIP_Cooper2014_rep2_ERR476761_ERR476757_trimV5_mm10.bam \
    --outdir cooper_k27me3
cp cooper_k27me3/NA_peaks.broadPeak ./ESC_serum_TKO_H3K27me3_ChIP_Cooper2014_broadPeak.bed

macs2 callpeak
    --broad \
    -g mm \
    --keep-dup all \
    --format BAM \
    -t ESC_serum_WT_H3K27me3_ChIP_Cooper2014_rep1_ERR476762_trimV5_mm10.bam ESC_serum_WT_H3K27me3_ChIP_Cooper2014_rep2_ERR476763_trimV5_mm10.bam \
    -c ESC_serum_WT_input_ChIP_Cooper2014_rep2_ERR476764_ERR476758_trimV5_mm10.bam \
    --outdir cooper_k27me3_WT
cp cooper_k27me3_WT/NA_peaks.broadPeak ./ESC_serum_WT_H3K27me3_ChIP_Cooper2014_broadPeak.bed

# SUZ12
macs2 callpeak
    --broad \
    -g mm \
    --keep-dup all \
    --format BAM \
    -t ESC_serum_TKO_Suz12_ChIP_Cooper2014_rep1_ERR476745_ERR476746_trimV5_mm10.bam ESC_serum_TKO_Suz12_ChIP_Cooper2014_rep2_ERR476750_ERR476772_trimV5_mm10.bam \
    -c ESC_serum_TKO_input_ChIP_Cooper2014_rep1_ERR476766_ERR476752_trimV5_mm10.bam  ESC_serum_TKO_input_ChIP_Cooper2014_rep2_ERR476761_ERR476757_trimV5_mm10.bam \
    --outdir cooper_Suz12
cp cooper_Suz12/NA_peaks.broadPeak ./ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak.bed

macs2 callpeak
    --broad \
    -g mm \
    --keep-dup all \
    --format BAM \
    -t ESC_serum_WT_Suz12_ChIP_Cooper2014_rep1_ERR476742_ERR476748_trimV5_mm10.bam ESC_serum_WT_input_ChIP_Cooper2014_rep2_ERR476764_ERR476758_trimV5_mm10.bam \
    -c ESC_serum_WT_input_ChIP_Cooper2014_rep2_ERR476764_ERR476758_trimV5_mm10.bam \
    --outdir cooper_Suz12_WT
cp cooper_Suz12_WT/NA_peaks.broadPeak ./ESC_serum_WT_Suz12_ChIP_Cooper2014_broadPeak.bed

# RING1B
macs2 callpeak
    --broad \
    -g mm \
    --keep-dup all \
    --format BAM \
    -t ESC_serum_TKO_Ring1B_ChIP_Cooper2014_rep1_ERR476756_ERR476770_trimV5_mm10.bam ESC_serum_TKO_Ring1B_ChIP_Cooper2014_rep2_ERR476759_ERR476765_trimV5_mm10.bam \
    -c ESC_serum_TKO_input_ChIP_Cooper2014_rep1_ERR476766_ERR476752_trimV5_mm10.bam  ESC_serum_TKO_input_ChIP_Cooper2014_rep2_ERR476761_ERR476757_trimV5_mm10.bam \
    --outdir cooper_Ring1b
cp cooper_Ring1b/NA_peaks.broadPeak ./ESC_serum_TKO_Ring1b_ChIP_Cooper2014_broadPeak.bed

macs2 callpeak
    --broad \
    -g mm \
    --keep-dup all \
    --format BAM \
    -t ESC_serum_WT_Ring1B_ChIP_Cooper2014_rep1_ERR476749_ERR476754_trimV5_mm10.bam ESC_serum_WT_Ring1B_ChIP_Cooper2014_rep2_ERR476755_ERR476753_trimV5_mm10.bam \
    -c ESC_serum_WT_input_ChIP_Cooper2014_rep2_ERR476764_ERR476758_trimV5_mm10.bam \
    --outdir cooper_Ring1b_WT
cp cooper_Ring1b_WT/NA_peaks.broadPeak ./ESC_serum_WT_Ring1b_ChIP_Cooper2014_broadPeak.bed



# stich together peaks within 10kb from files that start with PREFIX
PREFIX=$1
LIST=$(ls $PREFIX*bed)

cat $LIST | sort -k1,1 -k2,2n | grep -v "random" | grep -v "Un" | awk '{OFS=FS="\t"} { print $1, $2-10000, $3+10000}' > $PREFIX"_10kb_pad.bed"
bedtools merge -i $PREFIX"_10kb_pad.bed"  > $PREFIX"_10kb_pad_merge.bed"
awk '{OFS=FS="\t"}{ print $1, $2+10000, $3-10000}' $PREFIX"_10kb_pad_merge.bed" > $PREFIX"_peaks_10kb_stitch.bed"

# k-means clustering using Morpheus
awk '{OFS=FS="\t"}{ if ($4==1) print $0}' ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_data.txt > ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_k1.bed
awk '{OFS=FS="\t"}{ if ($4==2) print $0}' ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_data.txt > ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_k2.bed
awk '{OFS=FS="\t"}{ if ($4==3) print $0}' ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_data.txt > ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_k3.bed
awk '{OFS=FS="\t"}{ if ($4==4) print $0}' ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_data.txt > ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_k4.bed
awk '{OFS=FS="\t"}{ if ($4==5) print $0}' ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_data.txt > ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_k5.bed
awk '{OFS=FS="\t"}{ if ($4==6) print $0}' ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_data.txt > ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_k6.bed

# merge the k-means groups that correspond to SWRs
cat ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_k1.bed ESC_serum_TKO_Suz12_ChIP_Cooper2014_broadPeak_peaks_10kb_stitch_over1kb_k4.bed | sort -k1,1 -k2,2n | cut -f 1-3 > SWRs_Cooper2014.bed
cat H3K27me3_peaks_10kb_stitch_k6_k3.bed H3K27me3_peaks_10kb_stitch_k6_k5.bed | sort -k1,1 -k2,2n | cut -f1-3 > SWRs_revision_H3K27me3based.bed
cat H3K27me3_SEACR_peaks_merged_counts_RPKM_150RPKM_10866_calc_kmeans_k1.bed H3K27me3_SEACR_peaks_merged_counts_RPKM_150RPKM_10866_calc_kmeans_k4.bed | sort -k1,1 -k2,2n > SWRs_first_round.bed

# count number of peaks/domains
wc -l SWRs_* # subtract 1 for header
    1621 SWRs_Cooper2014.bed
    2374 SWRs_first_round.bed
    1302 SWRs_revision_H3K27me3based.bed

# I definitely had split domains in the original analysis. OK, redo everything.
cat H3K27me3_peaks_10kb_stitch_k6_k1.bed H3K27me3_peaks_10kb_stitch_k6_k4.bed H3K27me3_peaks_10kb_stitch_k6_k6.bed > CORs_revision_H3K27me3based.bed
cat H3K27me3_peaks_10kb_stitch_k6_k2.bed > ESRs_revision_H3K27me3based.bed

# forgot 10kb bins
cat CORs_revision_H3K27me3based.bed ESRs_revision_H3K27me3based.bed SWRs_revision_H3K27me3based_kmeans.bed CpG_islands_mm10_2.bed | sort -k1,1 -k2,2n > ROIs_revised.bed
