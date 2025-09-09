#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=250GB
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --job-name=correct2_final
#SBATCH --time=15:00:00
#SBATCH --partition=general
#SBATCH --account=a_nguyen_quan

source /etc/profile.d/modules.sh
module load miniconda3/4.12.0
source activate /home/s4716765/.conda/envs/scTAR_cellranger
cd /scratch/project_mnt/S0010/Prakrithi/lnc_revision/characterization_re/split/correct

awk 'BEGIN{OFS="\t"} {print $1, $10, $10+1, "donor_"$5"\t.\t"$4}' /QRISdata/Q4386/reference_anno_files/intropolis.v1.hg19_with_liftover_to_hg38.tsv | awk '!/NA/ {print}' > intropolis_donor_sites.bed
awk 'BEGIN{OFS="\t"} {print $1, $11-1, $11, "acceptor_"$6"\t.\t"$4}' /QRISdata/Q4386/reference_anno_files/intropolis.v1.hg19_with_liftover_to_hg38.tsv | awk '!/NA/ {print}' > intropolis_acceptor_sites.bed

# Input files
cuTARs_bed="/QRISdata/Q4386/all_cuTARs/all_cuTARs_refined_after_gencode47_clean_v5.bed"
cuTARs_5_prime_bed="cuTARs_5prime_new.bed"  # 5' end of cuTARs
cuTARs_3_prime_bed="cuTARs_3prime_new.bed"  # 3' end of cuTARs
TSS_bed="/QRISdata/Q4386/reference_anno_files/TSS_human_hg38_FANTOM_CAGE.bed"
intropolis_donors_bed="intropolis_donor_sites.bed"
intropolis_acceptors_bed="intropolis_acceptor_sites.bed"
polyA_sites_bed="/QRISdata/Q4386/reference_anno_files/PAS_atlas.clusters.3.0.GRCh38.GENCODE_42.bed"

# Output files and directory
output_summary="cuTAR_exon_classification_summary.txt"
output_summary_final="cuTAR_exon_classification_summary_final.txt"
output_dir="overlap_outputs_final"
mkdir -p "$output_dir"  # Create output directory if it doesn’t exist

# Temporary directory
tmp_dir=$(mktemp -d)

# Create header for the initial summary table
echo -e "cuTAR_ID\tClassification\tEvidence" > "$output_summary"

# Step 1: Full Transcript (TSS at 5' end AND PAS at 3' end)
echo "Finding full transcripts (TSS at 5' end and PAS at 3' end)..."
bedtools window -a "$cuTARs_5_prime_bed" -b "$TSS_bed" -l 1500 -r 100 -sm | \
    awk '{print $4}' | sort -u > "$tmp_dir/cutars_with_TSS.txt"
bedtools window -a "$cuTARs_3_prime_bed" -b "$polyA_sites_bed" -l 500 -r 1500 -sm | \
    awk '{print $4}' | sort -u > "$tmp_dir/cutars_with_PAS.txt"
comm -12 "$tmp_dir/cutars_with_TSS.txt" "$tmp_dir/cutars_with_PAS.txt" | \
    awk '{print $1"\tFull Transcript\tTSS,PAS"}' | sort -u > "$output_dir/full_transcript.bed"
cat "$output_dir/full_transcript.bed" >> "$output_summary"

# Step 2: First Exon (TSS at 5' end AND Donor at 3' end)
echo "Finding first exons (TSS at 5' end and Donor at 3' end)..."
bedtools window -a "$cuTARs_5_prime_bed" -b "$TSS_bed" -l 1500 -r 100 -sm > "$output_dir/cuTARs_5prime_with_TSS.bed"
bedtools window -a "$cuTARs_3_prime_bed" -b "$intropolis_donors_bed" -l 500 -r 500 -sm > "$output_dir/cuTARs_3prime_with_donor.bed"
awk 'NR==FNR {tss[$4]; next} $4 in tss' "$output_dir/cuTARs_5prime_with_TSS.bed" "$output_dir/cuTARs_3prime_with_donor.bed" | \
    sort -u -k4,4 | grep -v -f "$tmp_dir/cutars_with_TSS.txt" | \
    awk '{print $4"\tFirst Exon\tTSS,Donor Junction"}' | sort -u > "$output_dir/first_exon.bed"
cat "$output_dir/first_exon.bed" >> "$output_summary"

# Step 3: Internal Exon (Acceptor at 5' end AND Donor at 3' end)
echo "Finding internal exons (Acceptor at 5' end and Donor at 3' end)..."
bedtools window -a "$cuTARs_5_prime_bed" -b "$intropolis_acceptors_bed" -l 500 -r 500 -sm > "$output_dir/cuTARs_5prime_with_acceptor.bed"
awk 'NR==FNR {acceptor[$4]; next} $4 in acceptor' "$output_dir/cuTARs_5prime_with_acceptor.bed" "$output_dir/cuTARs_3prime_with_donor.bed" | \
    sort -u -k4,4 | grep -v -f "$tmp_dir/cutars_with_TSS.txt" | grep -v -f "$output_dir/first_exon.bed" | \
    awk '{print $4"\tInternal Exon\tAcceptor,Donor Junctions"}' | sort -u > "$output_dir/internal_exon.bed"
cat "$output_dir/internal_exon.bed" >> "$output_summary"

# Step 4: Last Exon (Acceptor at 5' end AND PAS at 3' end)
echo "Finding last exons (Acceptor at 5' end and PAS at 3' end)..."
bedtools window -a "$cuTARs_3_prime_bed" -b "$polyA_sites_bed" -l 500 -r 1500 -sm > "$output_dir/cuTARs_3prime_with_PAS.bed"
awk 'NR==FNR {acceptor[$4]; next} $4 in acceptor' "$output_dir/cuTARs_5prime_with_acceptor.bed" "$output_dir/cuTARs_3prime_with_PAS.bed" | \
    sort -u -k4,4 | grep -v -f "$tmp_dir/cutars_with_TSS.txt" | grep -v -f "$output_dir/first_exon.bed" | grep -v -f "$output_dir/internal_exon.bed" | \
    awk '{print $4"\tLast Exon\tAcceptor Junction,PAS"}' | sort -u > "$output_dir/last_exon.bed"
cat "$output_dir/last_exon.bed" >> "$output_summary"

# Step 5: No Evidence (No overlap with TSS, Donor, Acceptor, or PAS)
echo "Finding cuTARs with no evidence..."
awk '{print $4}' "$cuTARs_bed" | sort -u > "$tmp_dir/all_cuTARs.txt"
cat "$output_dir/full_transcript.bed" "$output_dir/first_exon.bed" "$output_dir/internal_exon.bed" "$output_dir/last_exon.bed" | \
    awk '{print $1}' | sort -u > "$tmp_dir/classified_cuTARs.txt"
comm -23 "$tmp_dir/all_cuTARs.txt" "$tmp_dir/classified_cuTARs.txt" | \
    awk '{print $1"\tNo Evidence\tNone"}' | sort -u > "$output_dir/no_evidence.bed"
cat "$output_dir/no_evidence.bed" >> "$output_summary"

# Move temporary files to output_dir
mv "$tmp_dir/cutars_with_TSS.txt" "$tmp_dir/cutars_with_PAS.txt" "$tmp_dir/all_cuTARs.txt" "$tmp_dir/classified_cuTARs.txt" "$output_dir/"

# Post-process the summary to ensure one classification per cuTAR
awk '
NR==1 {print "cuTAR_ID\tClassification\tEvidence"; next}  # Print header
{
    id = $1; class = $2; evid = $3;
    if (!(id in seen)) {  # First classification for this cuTAR
        if (class == "Full Transcript") {
            seen[id] = id "\t" class "\t" evid; priority[id] = 1;
        } else if (class == "First Exon") {
            seen[id] = id "\t" class "\t" evid; priority[id] = 2;
        } else if (class == "Last Exon") {
            seen[id] = id "\t" class "\t" evid; priority[id] = 3;
        } else if (class == "Internal Exon") {
            seen[id] = id "\t" class "\t" evid; priority[id] = 4;
        } else if (class == "No Evidence") {
            seen[id] = id "\t" class "\t" evid; priority[id] = 5;
        }
    }
}
END {
    for (id in seen) print seen[id];
}' "$output_summary" | sort -k1,1 > "$output_summary_final"

# Output the final result summary
echo "Initial summary saved to $output_summary"
echo "Final summary saved to $output_summary_final"
echo "All outputs saved in $output_dir"
echo "Temporary directory retained at $tmp_dir"
echo "Total cuTARs classified:"
wc -l "$output_summary_final"
echo "Classification counts:"
awk 'NR>1 {count[$2]++} END {for (type in count) print type "\t" count[type]}' "$output_summary_final" | sort
