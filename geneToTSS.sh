# JRA 2018
# Check the number of command line arguments

if [ $# -ne 1 ]; then
        script_name=$(basename $0)
        echo "Usage: $script_name input_bed"
        exit 1
fi

# REQUIREMENTS
# bedtools v.2.27.0
# bedops v.2.4.26
# sed GNU

INPUT_FILE=$1
FILE=$(basename "$INPUT_FILE")

# split by strand
awk '{OFS="\t";FS="\t"}{ if ($4 == "+") $3=$2+1; if ($4 == "+") print $0}' "$FILE" > "$FILE"_pos
awk '{OFS="\t";FS="\t"}{ if ($4 == "-") $2=$3-1; if ($4 == "-") print $0}' "$FILE" > "$FILE"_neg

# merge strands
cat "$FILE"_pos "$FILE"_neg > "$FILE"_unsort.bed

# sort by coordinate and prepend header
sort -V "$FILE"_unsort.bed > "$FILE"_sort.bed
touch header_for_geneToTSS
echo "#chr      start   end     strand	name" > header_for_geneToTSS
cat header_for_geneToTSS "$FILE"_sort.bed > ${FILE/%.*/_TSS.bed5}
echo "Finished creating bed5 file for $FILE"

rm "$FILE"_pos "$FILE"_neg "$FILE"_unsort.bed "$FILE"_sort.bed header_for_geneToTSS

COLUMN_COUNT=$(cat ${FILE/%.*/_TSS.bed5} | awk '{print NF}' | tail -n 1)
if [ $COLUMN_COUNT != 5 ]; then
	echo "Your file has $COLUMN_COUNT columns"
	echo "Please edit header appropriately"
fi

echo "$(basename $0) done!"
echo ""
echo ""
