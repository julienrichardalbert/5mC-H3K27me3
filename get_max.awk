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
} } END {print p}' SWRs_index_0U_0D_100w_50s_RING1B_Joshi2015.bed > SWRs_index_0U_0D_100w_50s_RING1B_Joshi2015_seed.bed
