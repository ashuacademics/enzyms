import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def extract_compound_name_from_filename(filename):
    # Extract the compound name from the input filename
    # Assuming the filename format is "processed_CompoundName.csv"
    parts = os.path.splitext(os.path.basename(filename))[0].split('_')
    if len(parts) == 2 and parts[0] == "processed":
        return parts[1]
    else:
        return None

def load_column_types(column_type_file):
    # Updated function to load sample columns based on the new file format with headers
    column_types = pd.read_csv(column_type_file)
    sample_columns = column_types[column_types['Sample'] == 'sample']['Enzyme'].tolist()
    return sample_columns

def plot_3d_mountain(input_files, column_type_file, output_directory):
    sample_columns = load_column_types(column_type_file)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for input_file in input_files:
        compound_name = extract_compound_name_from_filename(input_file)
        if compound_name is None:
            print(f"Invalid input file format: {input_file}")
            continue

        data = pd.read_csv(input_file)
        mz_values = data['mz'].copy()
        retention_time = data['rtime'].copy()

        # Filter data to include only relevant sample columns
        relevant_sample_columns = [col for col in data.columns if col in sample_columns]
        data = data[relevant_sample_columns]

        # Create a 3D plot
        fig = plt.figure(figsize=(15, 12))
        ax = fig.add_subplot(111, projection='3d')

        # Create meshgrid for X, Y, and determine Z
        x = np.arange(len(data.columns))
        y = np.arange(len(data.index))
        X, Y = np.meshgrid(x, y)
        Z = data.values

        # Plot 3D surface
        ax.plot_surface(X, Y, Z, cmap="coolwarm", edgecolor='k')

        # Labeling and ticks
        ax.set_xlabel('Samples', labelpad=50, fontsize=16)
        ax.set_ylabel('Retention Time', labelpad=20, fontsize=16)
        ax.set_zlabel('Peak Areas', labelpad=10, fontsize=16)

        # x-ticks
        num_xticks = min(10, len(data.columns))
        xticks = np.linspace(0, len(data.columns) - 1, num_xticks, dtype=int)
        ax.set_xticks(xticks)
        ax.set_xticklabels(data.columns[xticks], rotation=90)

        # y-ticks
        num_yticks = min(10, len(data.index))
        yticks = np.linspace(0, len(data.index) - 1, num_yticks, dtype=int)
        ax.set_yticks(yticks)
        ax.set_yticklabels(retention_time[yticks], rotation=90)

        # Save the plot to a PNG file
        file_base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file_path = os.path.join(output_directory, f"mountain_plot_{file_base_name}.png")
        plt.savefig(output_file_path)
        plt.close()

        print(f"Saved 3D mountain plot for file: {input_file}")

        # Prepare data of the top peaks with retention time
        top_peaks_data = []

        # Sort the peaks and select the top 100
        sorted_peak_values = np.sort(Z, axis=None)[::-1][:100]
        sorted_peak_indices = np.argwhere(np.isin(Z, sorted_peak_values))

        for peak_idx in sorted_peak_indices:
            row_idx, col_idx = peak_idx[0], peak_idx[1]  # Corrected the typo here
            sample_name = data.columns[col_idx]
            mz_value = mz_values.iloc[row_idx]
            rt_value = retention_time.iloc[row_idx]
            peak_value = Z[row_idx, col_idx]

            top_peaks_data.append({
            'Sample': sample_name,
            'M/Z': mz_value,
            'Retention Time (sec)': rt_value,
            'Value': peak_value
        })


        # Define the output filename for the highest peaks data
        peaks_output_file_name = f"highest_peaks_{file_base_name}.csv"

        # Save the highest peaks to a CSV file
        peaks_df = pd.DataFrame(top_peaks_data)
        peaks_output_file_path = os.path.join(output_directory, peaks_output_file_name)
        peaks_df.to_csv(peaks_output_file_path, index=False)
        print(f"Saved CSV file of the highest peaks for file: {input_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot 3D mountain plots from CSV files.')
    parser.add_argument('--input_files', type=str, nargs='+', help='Paths to the input CSV files')
    parser.add_argument('--column_type_file', type=str, help='Path to the file with column types')
    parser.add_argument('--output_directory', type=str, help='Directory to save the 3D mountain plots and peak data CSVs')

    args = parser.parse_args()
    plot_3d_mountain(args.input_files, args.column_type_file, args.output_directory)
