echo '++++++++++++++++++++++++++++++++++++++++++++++++'
echo '+++++         Enriched Motifs              +++++'
echo '++++++++++++++++++++++++++++++++++++++++++++++++'


bed_dir="biomarker/atac_bio_up/bed_files/Chen-2019-small"
output_dir="biomarker/atac_bio_up/enriched_motifs_other/Chen-2019-small"
back_dir="biomarker/atac_bio_up/back_files_other/Chen-2019-small"

for bed_file in "$bed_dir"/*.bed; do
    filename=$(basename "$bed_file" .bed)

    celltype=$(echo "$filename" | cut -d'_' -f1)
    tool=$(echo "$filename" | cut -d'_' -f2-)

    echo Processing $celltype with tool $tool

    gimme background "$back_dir/${filename}.bed" gc -i $bed_file -f BED -n 1000 -g mm10.fa
    gimme motifs $bed_file "$output_dir/${filename}/" --known -g mm10.fa -b "$back_dir/${filename}.bed" --nogc
done