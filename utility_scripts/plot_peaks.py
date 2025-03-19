#/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
from pyopenms import MSExperiment, MzMLFile
import pandas as pd
import os

'''
This script plots the intensity profile of a specified product from multiple mzML files.
The script requires the following arguments:
1. -i, --input: Paths to mzML files.
2. -p, --product: m/z of the product.
3. -t, --tolerance: m/z tolerance for peak matching.
4. -o, --output: Output path for the PNG plot.
5. -s, --sample_file: Path to the sample mapping file.
'''

def read_sample_mapping(sample_file):
    sample_df = pd.read_csv(sample_file)
    return sample_df.set_index('Sample')['Enzyme'].to_dict()

def plot_product_peaks(mzml_paths, mz_product, tolerance, output_path, sample_mapping):
    # Initialize data structures for plotting
    rt_products = []
    intensity_products = []
    sample_names = []
    colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan', 'magenta', 'yellow', 'black', 'brown']  # Extended color list

    # Function to extract data
    def extract_data(exp, mz_value):
        rt_list, intensity_list = [], []  # Initialize the lists properly
        for spectrum in exp:
            if spectrum.getMSLevel() == 1:  # Only consider MS1 spectra
                rt = spectrum.getRT()
                mz_array, intensity_array = spectrum.get_peaks()
                for mz, intensity in zip(mz_array, intensity_array):
                    if abs(mz - mz_value) <= tolerance:
                        rt_list.append(rt)
                        intensity_list.append(intensity)
        return rt_list, intensity_list

    # Extract data for each file and prepare sample names
    for mzml_path in mzml_paths:
        exp = MSExperiment()
        MzMLFile().load(mzml_path, exp)
        rt_product, intensity_product = extract_data(exp, mz_product)
        rt_products.append(rt_product)
        intensity_products.append(intensity_product)

        # Extract and format the sample name from the sample mapping file
        file_name = os.path.basename(mzml_path)
        sample_name_key = file_name.replace('.mzML', '').split('_')[1]  # Select the prefix
        sample_name = sample_mapping.get(sample_name_key, sample_name_key)
        sample_names.append(sample_name)

    # Determine the number of rows needed for two columns
    num_samples = len(mzml_paths)
    num_columns = 1
    num_rows = (num_samples + num_columns - 1) // num_columns

    # Plotting
    fig, axs = plt.subplots(num_rows, num_columns, figsize=(10, 4 * num_rows), sharex=True)

    # Flatten the axs array for easy iteration
    axs = axs.flatten()

    # Plot for each sample
    for i, (rt_product, intensity_product, sample_name) in enumerate(zip(rt_products, intensity_products, sample_names)):
        if rt_product and intensity_product:
            label = f'{sample_name}'
            color = colors[i % len(colors)]
            axs[i].plot(rt_product, intensity_product, '-', linewidth=2, alpha=0.7, color=color)
            axs[i].scatter(rt_product, intensity_product, color=color, s=10, alpha=0.7)
            axs[i].fill_between(rt_product, intensity_product, color=color, alpha=0.3)
            axs[i].set_ylabel('Peak Intensity', fontsize=28)
            
            axs[i].set_xticks(range(60, 71, 2))
            axs[i].tick_params(axis='both', which='major', labelsize=24)
            legend = axs[i].legend([label], fontsize=24, frameon=False)
            for handle in legend.legend_handles:
                handle.set_visible(False)
            axs[i].set_ylim(0, 200000)
            axs[i].set_xlim(60, 70)
        else:
            axs[i].remove()  # Remove empty subplots

    # Finalizing the plot
    axs[-1].set_xlabel('Retention Time (sec)', fontsize=18)
    plt.tight_layout()
    plt.savefig(output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot intensity profile for specified product from multiple mzML files.')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Paths to mzML files.')
    parser.add_argument('-p', '--product', type=float, required=True, help='m/z of product.')
    parser.add_argument('-t', '--tolerance', type=float, default=0.1, help='m/z tolerance for peak matching.')
    parser.add_argument('-o', '--output', required=True, help='Output path for the PNG plot.')
    parser.add_argument('-s', '--sample_file', required=True, help='Path to the sample mapping file.')

    args = parser.parse_args()

    sample_mapping = read_sample_mapping(args.sample_file)
    plot_product_peaks(args.input, args.product, args.tolerance, args.output, sample_mapping)
