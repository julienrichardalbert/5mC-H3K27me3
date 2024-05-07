#!/bin/bash
INPUT=$1

#head -n 3 ESRs_kmeans6_4569_1000w_50s_data_100bp.bed
#chr1	17096413	17097413	1	19.767723089547665	6.280056356957005	17.282458104150702	
#18.63474196223893

awk '{
if($4 == val) {
if ( $5+$6 > max ){
p=$0;
max = $5+$6;
}
} else {
print p;
val=$4;
max=$5+$6;
p=$0;
} } END {print p}' $INPUT > ${INPUT//.bed/_RING1B_Joshi2015_seed.bed}


awk '{
if($4 == val) {
if ( $7+$8 > max ){
p=$0;
max = $7+$8;
}
} else {
print p;
val=$4;
max=$7+$8;
p=$0;
} } END {print p}' $INPUT > ${INPUT//.bed/_ESC_H3K27me3_CnT_seed.bed}
