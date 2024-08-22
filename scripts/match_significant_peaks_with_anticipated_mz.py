import pandas as pd
import os
import argparse

def load_reference(reference_path):
    return pd.read_csv(reference_path)

def extract_sample_name_from_filename(filename):
    # Extracts the sample name from the filename, assuming the format "highest_peaks_processed_SampleName.csv"
    return os.path.splitext(os.path.basename(filename))[0].split('_')[-1]

def find_matching_mz_values(sample_path, sample_name, reference_data, tolerance=0.01):
    sample_data = pd.read_csv(sample_path)

    matched_rows = []

    for _, sample_row in sample_data.iterrows():
        input_sample_name = sample_row['Sample']  # Get the sample name from the input file
        mz_value = sample_row['M/Z']
        retention_time = sample_row['Retention Time (sec)']
        value = sample_row['Value']

        for _, reference_row in reference_data.iterrows():
            substrate = reference_row['Substrate']
            simulated_mz = reference_row['Simulated-m/z']
            simulated_formula = reference_row['SimulatedFormula']

            difference = abs(mz_value - simulated_mz)

            if difference <= tolerance and substrate == sample_name:
                matched_rows.append({
                    'Sample': input_sample_name,  # Include the original sample name
                    'M/Z': mz_value,
                    'Substrate': substrate,
                    'Simulated-m/z': simulated_mz,
                    'SimulatedFormula': simulated_formula,
                    'Retention Time (sec)': retention_time,
                    'Value': value
                })
                break

    return pd.DataFrame(matched_rows)



def process_files(input_dir, reference_file, output_dir, tolerance):
    reference_data = load_reference(reference_file)

    for filename in os.listdir(input_dir):
        if filename.startswith("highest_peaks_processed_") and filename.endswith(".csv"):
            sample_name = extract_sample_name_from_filename(filename)
            data_path = os.path.join(input_dir, filename)

            matched_data = find_matching_mz_values(data_path, sample_name, reference_data, tolerance)

            if not matched_data.empty:
                output_filename = f"output_{sample_name}.csv"
                output_path = os.path.join(output_dir, output_filename)
                matched_data.to_csv(output_path, index=False)
                print(f"Written matches for {sample_name} to {output_path}")
            else:
                print(f"No matches found for {sample_name}.")
        else:
            print(f"Skipping {filename} as it does not match the expected naming convention.")

def main():
    parser = argparse.ArgumentParser(description='Filter matching M/Z values.')
    parser.add_argument('--input_dir', type=str, help='Directory containing the input CSV files')
    parser.add_argument('--reference_file', type=str, help='Path to the reference CSV file')
    parser.add_argument('--output_dir', type=str, help='Directory where output CSV files will be saved')
    parser.add_argument('--tolerance', type=float, default=0.01, help='Tolerance for matching (default: 0.01)')
    
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    process_files(args.input_dir, args.reference_file, args.output_dir, args.tolerance)

if __name__ == "__main__":
    main()
