#/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
from pyopenms import MSExperiment, MzMLFile

'''
Plot intensity profile for specified substrate and product from mzML file and compare with empty vector control.

Usage:
python plot_substrate_product_peak_with_evc.py -i <main_mzml_file> -c <control_mzml_file> -s <substrate_mz> -p <product_mz> -t <tolerance> -o <output_png_file>

'''

def plot_filtered_peaks_with_fill(mzml_path, control_mzml_path, mz_substrate, mz_product, tolerance, output_path):
    # Load mzML file for the main sample
    exp = MSExperiment()
    MzMLFile().load(mzml_path, exp)
    
    # Load mzML file for the control sample
    control_exp = MSExperiment()
    MzMLFile().load(control_mzml_path, control_exp)

    # Prepare data structures for plotting
    rt_substrate, intensity_substrate, rt_product, intensity_product = [], [], [], []
    rt_substrate_control, intensity_substrate_control, rt_product_control, intensity_product_control = [], [], [], []

    # Function to extract data
    def extract_data(exp, rt_list, intensity_list, mz_value):
        for spectrum in exp:
            if spectrum.getMSLevel() == 1:  # Only consider MS1 spectra
                rt = spectrum.getRT()
                mz_array, intensity_array = spectrum.get_peaks()
                for mz, intensity in zip(mz_array, intensity_array):
                    if abs(mz - mz_value) <= tolerance:
                        rt_list.append(rt)
                        intensity_list.append(intensity)

    # Extract data for main sample
    extract_data(exp, rt_substrate, intensity_substrate, mz_substrate)
    extract_data(exp, rt_product, intensity_product, mz_product)

    # Extract data for control sample
    extract_data(control_exp, rt_substrate_control, intensity_substrate_control, mz_substrate)
    extract_data(control_exp, rt_product_control, intensity_product_control, mz_product)

    # Plotting
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)

    # Plot for main sample with points
    if rt_substrate and intensity_substrate:
        ax.plot(rt_substrate, intensity_substrate, 'b-', label=f'Sample Substrate', linewidth=2, alpha=0.7)
        #ax.scatter(rt_substrate, intensity_substrate, color='blue', s=10, alpha=0.7)  # Plot points only for main sample
        ax.fill_between(rt_substrate, intensity_substrate, color='blue', alpha=0.3)
    if rt_product and intensity_product:
        ax.plot(rt_product, intensity_product, 'r-', label=f'Sample Product', linewidth=2, alpha=0.7)
        #ax.scatter(rt_product, intensity_product, color='darkred', s=10, alpha=0.7)  # Plot points only for main sample
        ax.fill_between(rt_product, intensity_product, color='red', alpha=0.3)

    # Plot for control sample without points
    if rt_substrate_control and intensity_substrate_control:
        ax.plot(rt_substrate_control, intensity_substrate_control, 'b--', label=f'EVC Substrate', linewidth=2, alpha=0.7)
        #ax.fill_between(rt_substrate_control, intensity_substrate_control, color='lightblue', alpha=0.3)
    if rt_product_control and intensity_product_control:
        ax.plot(rt_product_control, intensity_product_control, 'r--', label=f'EVC Product', linewidth=2, alpha=0.7)
        #ax.fill_between(rt_product_control, intensity_product_control, color='pink', alpha=0.3)

    # Finalizing the plot
    ax.set_xlabel('Retention Time (sec)', fontsize=28)
    ax.set_ylabel('Peak Intensity', fontsize=28)
    ax.legend(fontsize=20, loc='upper left', frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=24)
    #plt.title('Peak Intensity Over Retention Time', fontsize=18)
    plt.xlim(60, 75)
    plt.tight_layout()
    plt.savefig(output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot intensity profile with filled peaks for specified substrate and product from mzML file.')
    parser.add_argument('-i', '--input', required=True, help='Path to the main mzML file.')
    parser.add_argument('-c', '--control', required=True, help='Path to the control mzML file.')
    parser.add_argument('-s', '--substrate', type=float, required=True, help='m/z of substrate.')
    parser.add_argument('-p', '--product', type=float, required=True, help='m/z of product.')
    parser.add_argument('-t', '--tolerance', type=float, default=0.1, help='m/z tolerance for peak matching.')
    parser.add_argument('-o', '--output', required=True, help='Output path for the PNG plot.')

    args = parser.parse_args()

    plot_filtered_peaks_with_fill(args.input, args.control, args.substrate, args.product, args.tolerance, args.output)
