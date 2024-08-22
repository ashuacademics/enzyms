import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import glob
import numpy as np
from matplotlib.ticker import FuncFormatter

def detect_isotope_pairs(df, delta=1.00335, tolerance=0.005):
    unique_mz_values = df['M/Z'].unique()
    pairs = []
    for mz in unique_mz_values:
        matching_mz = unique_mz_values[np.abs(unique_mz_values - (mz + delta)) <= tolerance]
        for match in matching_mz:
            pairs.append((mz, match))
    return pairs

def read_substrate_mz(file_path):
    substrate_df = pd.read_csv(file_path)
    return substrate_df.set_index('Substrate')['Substrate-m/z'].to_dict()

def read_sample_mapping(sample_file):
    sample_df = pd.read_csv(sample_file)
    return sample_df.set_index('Sample')['Enzyme'].to_dict()

def remove_isotopes_of_substrate(df, substrate_mz_dict, delta=1.00335, tolerance=0.005):
    isotopes_to_remove = set()
    for substrate, mz in substrate_mz_dict.items():
        isotope_mz = mz + delta
        isotopes_to_remove.update(df['M/Z'][np.abs(df['M/Z'] - isotope_mz) <= tolerance])
    return df[~df['M/Z'].isin(isotopes_to_remove)]

def plot_scatter(csv_file, output_dir, substrate_mz_dict, sample_mapping):
    df = pd.read_csv(csv_file)
    
    # Print the columns to debug the issue
    print("Columns in the CSV file:", df.columns)

    # Detect isotope pairs
    isotope_pairs = detect_isotope_pairs(df)
    isotopes_to_remove = set([mz2 for _, mz2 in isotope_pairs])

    # Filter out isotope pairs
    df_filtered = df[~df['M/Z'].isin(isotopes_to_remove)]
    
    # Remove isotopes of substrate
    df_filtered = remove_isotopes_of_substrate(df_filtered, substrate_mz_dict)

    # Normalize retention times to use as point sizes
    min_size = 50
    max_size = 200
    retention_time = df_filtered['Retention Time (sec)']
    #sizes = min_size + (retention_time - retention_time.min()) / (retention_time.max() - retention_time.min()) * (max_size - min_size)

    plt.figure(figsize=(12, 10))

    unique_mz_values = df_filtered['M/Z'].unique()
    unique_mz_values.sort()
    mz_to_x_position = {mz: idx for idx, mz in enumerate(unique_mz_values)}

    # Define marker styles and colors
    markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h', 'H', '+', 'x', 'd', '|', '_']
    colormap = plt.colormaps.get_cmap('tab20')

    for idx, sample in enumerate(df_filtered['Sample'].unique()):
        subset = df_filtered[df_filtered['Sample'] == sample]
        sample_suffix = sample_mapping.get(sample, sample)
        x_positions = [mz_to_x_position[mz] for mz in subset['M/Z']]
        plt.scatter(x_positions, subset['Value'], label=sample_suffix, alpha=1, s=150, #s=sizes[subset.index], 
                    marker=markers[idx % len(markers)], color=colormap(idx % 20))

    plt.xlabel('Mass-to-Charge Ratio (M/Z)', fontsize=24)
    plt.ylabel('Peak Intensity', fontsize=24)
    plt.xticks(range(len(unique_mz_values)), [str(x) for x in unique_mz_values], rotation=90)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    # legend font size
    #plt.legend(prop={'size': 14})


    ax = plt.gca()
    ax.yaxis.set_major_formatter(FuncFormatter(lambda val, pos: f'{int(val / 1e6)}M' if val >= 1e6 else f'{int(val / 1e3)}K' if val >= 1e3 else f'{int(val)}'))

    input_base_name = os.path.splitext(os.path.basename(csv_file))[0]
    output_file = os.path.join(output_dir, f'scatter_plot_{input_base_name}.png')
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Generate scatter plots from CSV files in a directory.')
    parser.add_argument('--input_dir', help='Path to the directory containing input CSV files')
    parser.add_argument('--output_dir', help='Path to the output directory')
    parser.add_argument('--substrate_mz_file', help='Path to the substrate m/z file')
    parser.add_argument('--sample_file', help='Path to the sample file')
    args = parser.parse_args()

    substrate_mz_dict = read_substrate_mz(args.substrate_mz_file)
    sample_mapping = read_sample_mapping(args.sample_file)

    input_files = glob.glob(os.path.join(args.input_dir, 'output_*.csv'))
    for input_file in input_files:
        plot_scatter(input_file, args.output_dir, substrate_mz_dict, sample_mapping)

if __name__ == "__main__":
    main()
