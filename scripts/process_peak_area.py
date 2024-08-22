import os
import argparse
import pandas as pd

def extract_compound_name(filename):
    return filename.split('_')[0]

def load_compound_mz(file):
    data = pd.read_csv(file)
    compound_mz_map = pd.Series(data['Substrate-m/z'].values, index=data['Substrate']).to_dict()
    return compound_mz_map

def load_column_types(column_type_file):
    column_types = pd.read_csv(column_type_file)
    control_columns = column_types[column_types['Sample'] == 'evc']['Enzyme'].tolist()
    sample_columns = column_types[column_types['Sample'] == 'sample']['Enzyme'].tolist()
    return control_columns, sample_columns

def process_csv(file, compound_name, compound_mz_map, control_columns, sample_columns, output_directory):
    if compound_name in compound_mz_map:
        compound_mz = compound_mz_map[compound_name]
        data = pd.read_csv(file, sep='\t')

        relevant_columns = control_columns + sample_columns

        for col in relevant_columns:
            data[col] = pd.to_numeric(data[col], errors='coerce')

        # Remove rows matching the m/z of the compound with a tolerance of 0.01
        tolerance = 0.01
        removed_rows = data[abs(data['mz'] - compound_mz) <= tolerance]
        data = data[abs(data['mz'] - compound_mz) > tolerance]

        # Remove the mean and std calculation for control columns when there's only one control column
        if len(control_columns) == 1:
            control_col = control_columns[0]
            for sample_col in sample_columns:
                data[sample_col] = (data[sample_col] - data[control_col]).clip(lower=0)
        else:
            control_mean = data[control_columns].mean(axis=1)
            control_std = data[control_columns].std(axis=1)

            for sample_col in sample_columns:
                data[sample_col] = (data[sample_col] - (control_mean + 3 * control_std)).clip(lower=0)

        # Filter the data as described
        filtered_data = data[
            (data['mz'] > (compound_mz - 100)) & 
            (data['mz'] < (compound_mz + 100)) & 
            ~((data['mz'] >= (compound_mz - 0.05)) & (data['mz'] <= (compound_mz + 0.05)))
        ]

        # Add the removed rows back to the processed data
        data = pd.concat([filtered_data, removed_rows])

        # Save the processed and filtered data to a new CSV file in the specified output directory
        output_file = os.path.join(output_directory, f"processed_{compound_name}.csv")
        filtered_data.to_csv(output_file, index=False)
        print(f"Processed data saved to {output_file}")

        # Create an additional output file with transposed data containing peak areas
        transposed_data = data[sample_columns + ['mz']].set_index('mz').T
        # Rename the index to 'Samples' for the transposed data
        transposed_data.index.name = 'Samples'
        transposed_output_file = os.path.join(output_directory, f"transposed_{compound_name}.csv")
        transposed_data.reset_index().to_csv(transposed_output_file, index=False)
        print(f"Transposed data saved to {transposed_output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process CSV files with argparse")
    parser.add_argument('--csv_directory', required=True, help='Directory containing the CSV files')
    parser.add_argument('--compound_mz_file', required=True, help='Path to the compound m/z file')
    parser.add_argument('--column_type_file', required=True, help='Path to the file with column types')
    parser.add_argument('--output_directory', required=True, help='Directory where processed files will be saved')

    args = parser.parse_args()

    compound_mz_map = load_compound_mz(args.compound_mz_file)
    control_columns, sample_columns = load_column_types(args.column_type_file)

    for file in os.listdir(args.csv_directory):
        if file.endswith('.csv'):
            file_path = os.path.join(args.csv_directory, file)
            compound_name = extract_compound_name(file)
            process_csv(file_path, compound_name, compound_mz_map, control_columns, sample_columns, args.output_directory)
