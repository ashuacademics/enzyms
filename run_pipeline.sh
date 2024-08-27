#!/bin/bash
# File              : run_pipeline.sh
# Author            : Ashutosh Kumar <kums@zhaw.ch>
# Date              : 27.08.2024
# Last Modified Date: 27.08.2024
# Last Modified By  : Ashutosh Kumar <kums@zhaw.ch>

# Ensure all necessary arguments are provided
if [ $# -ne 10 ]; then
    echo "Usage: $0 --smi_file <smi_file> --variations_file <variations_param_file> --params_file <asari_parameters_file> --samples_file <list_of_samples_file> --mzml_dir <mzML_directory>"
    exit 1
fi

# Parsing named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --smi_file) smi_file="$2"; shift ;;
        --variations_file) variations_param_file="$2"; shift ;;
        --params_file) asari_params_file="$2"; shift ;;
        --samples_file) list_of_samples_file="$2"; shift ;;
        --mzml_dir) mzml_dir="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Derive the compound name from the .smi file
compound_name=$(basename "$smi_file" .smi)

# Detect and align features using Asari
asari process --mode pos --input "${mzml_dir}" --reference "${mzml_dir}/${compound_name}_EVC1.mzML" -o ${compound_name} --parameter "$asari_params_file"

# Copy peak area matrix from Asari results
cp ${compound_name}_asari_project*/export/full_Feature_table.tsv ${compound_name}_areas.csv

# Calculate m/z of substrate
python3 scripts/calculate_substrate_mz.py --input_file ${compound_name}.smi --output_file substrate_mz.csv --parameters_file "$variations_param_file"

# Calculate m/z of anticipated products
python3 scripts/calculate_anticipated_products_mz.py --input_file ${compound_name}.smi --output_file Anticipated_products_mz.txt --parameters_file "$variations_param_file"

# Adjust peak areas
python3 scripts/process_peak_area.py --csv_directory ./ --compound_mz_file substrate_mz.csv --column_type_file "$list_of_samples_file" --output_directory ./

# Get top 100 peaks
python3 scripts/get_significant_peaks.py --input_files ./processed_${compound_name}.csv --column_type_file "$list_of_samples_file" --output_directory ./

# Match simulated products file
python3 scripts/match_significant_peaks_with_anticipated_mz.py --input_dir ./ --reference_file Anticipated_products_mz.txt --output_dir ./ --tolerance 0.01

# Plot pipeline scatterplot
python3 scripts/plot_scatterplot_matched_peaks.py --input_dir ./ --output_dir ./ --substrate_mz_file substrate_mz.csv --sample_file "$list_of_samples_file"

# Ensure the output directory exists
mkdir -p output

# Move all CSV and PNG files to the output directory
find . -maxdepth 1 -type f \( -name "*.csv" -o -name "*.png" \) -exec mv {} output/ \;

# Ensure the files in the output directory have the correct permissions to be removed later
chmod -R 777 output

echo "EnzyMS Pipeline execution completed for compound: $compound_name"
