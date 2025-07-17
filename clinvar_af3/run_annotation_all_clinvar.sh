#!/bin/bash
input_dir="af3_clinvar_ref"
reference_dir="../reference"
structure_informed_rvas_dir="../"
protein_list=$input_dir/protein_list.txt

for input in $input_dir/annotation_files_by_chr/*.tsv; do
    echo "Processing $input"
    
    chr_base=$(basename "$input" .tsv)
    output_dir="result/${chr_base}"
    mkdir -p "$output_dir"

    chr_num=$(basename "$input" | sed -E 's/^.*chr([0-9]+).*$/\1/')

    rvas_input="$input_dir/uniprot_clinvar_rvas/clinvar_grch38_annotated_2_chr${chr_num}_rvas.tsv.gz"

    if [[ -f "$rvas_input" ]]; then
        python "$structure_informed_rvas_dir/run.py" \
            --rvas-data-to-map "$rvas_input" \
            --reference-dir "$reference_dir" \
            --which-proteins "$protein_list" \
            --annotation_file "$input" \
            --results-dir "$output_dir/"
    else
        echo "Skipping chromosome ${chr_num} — missing $rvas_input"
    fi
done
