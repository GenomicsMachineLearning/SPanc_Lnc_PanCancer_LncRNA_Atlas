#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=250GB
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --job-name=overlap_grok_all
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_nguyen_quan

source /etc/profile.d/modules.sh
module load miniconda3/4.12.0
source activate /home/s4716765/.conda/envs/scTAR_cellranger

# List of sample .bed files
#samples=("sample1.bed" "sample2.bed" "sample3.bed")  # Update with actual filenames
#samples=($(ls /QRISdata/Q4386/lnc_revision/fig1/data/*.bed))

samples="/QRISdata/Q4386/all_cuTARs/all_cuTARs_refined_after_gencode47_clean_v5.bed"
# Define paths

# Directories and reference files (update paths as needed)
SAMPLE_DIR="/QRISdata/Q4386/lnc_revision/fig1/data/"
REF_FANTOM="/QRISdata/Q4386/lncrna_public_data/FANTOM_lncRNAs_hg38.bed"
REF_LNCBOOK="/QRISdata/Q4386/lncrna_public_data/LncBook_Version2.0_onlylnc.bed"
REF_NONCODE="/QRISdata/Q4386/reference_anno_files/NONCODEv5_hg38.lncAndGene.bed"
REF_GENCODE="/QRISdata/Q4386/reference_anno_files/gencode_v47_lncRNAonly.bed"
REF_LNCPEDIA="/QRISdata/Q4386/lncrna_public_data/lncipedia_5_2_hc_hg38.bed"

# Output file
OUTPUT="overlap_summary_grok_all_cuTARs.tsv" #overlap_summary_grok.tsv

TMP_DIR="./debug_all_cuTARs" # "./debug_temp2"  # Persistent debug directory

# Clean and create debug dir
rm -rf "$TMP_DIR" && mkdir -p "$TMP_DIR"

# Header
echo -e "Sample\tNovel\tFANTOM\tLncBook\tNONCODE\tGENCODE\tLncPedia\tIn_Two\tIn_Three\tIn_Four\tIn_Five" > "$OUTPUT"

# Process each sample
for SAMPLE in "$samples"; do # for only cutar bed
#for SAMPLE in "$SAMPLE_DIR"/*.bed; do #for sample-wise
    if [ ! -f "$SAMPLE" ]; then
        echo "No BED files found in $SAMPLE_DIR"
        exit 1
    fi

    SAMPLE_NAME=$(basename "$SAMPLE")
    TOTAL=$(wc -l < "$SAMPLE")
    echo "Processing $SAMPLE_NAME with $TOTAL regions"

    # Verify sample BED (ensure strand column exists)
    echo "Sample head (check strand):"
    head -n 5 "$SAMPLE"

    # Generate strand-specific overlaps
    REFS=("$REF_FANTOM" "$REF_LNCBOOK" "$REF_NONCODE" "$REF_GENCODE" "$REF_LNCPEDIA")
    DB_NAMES=("fantom" "lncbook" "noncode" "gencode" "lncpédia")
    for i in "${!REFS[@]}"; do
        DB="${DB_NAMES[$i]}"
        bedtools intersect -a "$SAMPLE" -b "${REFS[$i]}" -u -s > "$TMP_DIR/${SAMPLE_NAME}_${DB}_overlaps.bed"
        COUNT=$(wc -l < "$TMP_DIR/${SAMPLE_NAME}_${DB}_overlaps.bed")
        echo "Strand-specific overlaps with $DB: $COUNT"
        [ "$COUNT" -gt 0 ] && head -n 5 "$TMP_DIR/${SAMPLE_NAME}_${DB}_overlaps.bed"
    done

    # Add unique IDs
    awk '{print $0"\t"$1":"$2"-"$3}' "$SAMPLE" > "$TMP_DIR/${SAMPLE_NAME}_with_id.bed"
    echo "Sample with ID head:"
    head -n 5 "$TMP_DIR/${SAMPLE_NAME}_with_id.bed"

    # Track overlaps per region
    > "$TMP_DIR/${SAMPLE_NAME}_overlaps_with_id.txt"  # Clear file
    for DB in "${DB_NAMES[@]}"; do
        if [ -s "$TMP_DIR/${SAMPLE_NAME}_${DB}_overlaps.bed" ]; then
            bedtools intersect -a "$TMP_DIR/${SAMPLE_NAME}_with_id.bed" -b "$TMP_DIR/${SAMPLE_NAME}_${DB}_overlaps.bed" -u -s | \
            awk -v db="$DB" '{print $NF"\t"db}' >> "$TMP_DIR/${SAMPLE_NAME}_overlaps_with_id.txt"
        fi
    done
    if [ -s "$TMP_DIR/${SAMPLE_NAME}_overlaps_with_id.txt" ]; then
        echo "Overlaps with ID: $(wc -l < "$TMP_DIR/${SAMPLE_NAME}_overlaps_with_id.txt")"
        sort "$TMP_DIR/${SAMPLE_NAME}_overlaps_with_id.txt" | uniq > "$TMP_DIR/${SAMPLE_NAME}_unique_overlaps.txt"
        echo "Unique overlaps head:"
        head -n 5 "$TMP_DIR/${SAMPLE_NAME}_unique_overlaps.txt"
    else
        echo "No overlaps detected in $SAMPLE_NAME"
    fi

    # Count overlaps per region
    awk '
    BEGIN {
        novel = 0;
        single["fantom"] = 0; single["lncbook"] = 0; single["noncode"] = 0;
        single["gencode"] = 0; single["lncpédia"] = 0;
        in_two = 0; in_three = 0; in_four = 0; in_five = 0;
    }
    NR==FNR {
        overlaps[$1] = overlaps[$1] ? overlaps[$1]","$2 : $2;
        next
    }
    {
        id = $NF;
        if (!overlaps[id]) { novel++ }
        else {
            split(overlaps[id], dbs, ",");
            count = length(dbs);
            if (count == 1) { single[dbs[1]]++ }
            else if (count == 2) { in_two++ }
            else if (count == 3) { in_three++ }
            else if (count == 4) { in_four++ }
            else if (count == 5) { in_five++ }
        }
    }
    END {
        printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            sample, novel,
            single["fantom"], single["lncbook"], single["noncode"],
            single["gencode"], single["lncpédia"],
            in_two, in_three, in_four, in_five
    }' sample="$SAMPLE_NAME" "$TMP_DIR/${SAMPLE_NAME}_unique_overlaps.txt" "$TMP_DIR/${SAMPLE_NAME}_with_id.bed" >> "$OUTPUT"

done

echo "Results written to $OUTPUT"
echo "Debug files saved in $TMP_DIR"
