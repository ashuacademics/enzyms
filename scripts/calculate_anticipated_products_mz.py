#!/usr/bin/env python3

import yaml
import argparse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from pyopenms import EmpiricalFormula
import re
from itertools import product
import adducts  

'''
This script calculates the m/z of the anticipated products of a substrate based on its SMILES data. 
The script reads a .smi file containing the SMILES data of the substrates and their names. 
It then calculates the m/z of the anticipated products of the substrates by simulating variations in 
the number of Hydrogen, Carbon, Oxygen, Nitrogen, Fluorine, and Chlorine atoms in the substrate. 
The script writes the results to an output file.

Usage:

python calculate_anticipated_products_mz.py --input_file input.smi --output_file output.csv --parameters_file variations.param

Arguments:
--input_file: Path to the input .smi file
--output_file: Path to save the output results
--parameters_file: Path to the variations.param file

The variations.param file should include:
adduct: The adduct ion type
mode: The ionization mode (positive or negative)
H_variation: Variation range for Hydrogen (list of two integers)
C_variation: Variation range for Carbon (list of two integers)
O_variation: Variation range for Oxygen (list of two integers)
N_variation: Variation range for Nitrogen (list of two integers)
F_variation: Variation range for Fluorine (list of two integers)
Cl_variation: Variation range for Chlorine (list of two integers)
'''

def calculateMass(formula: str, adduct_mass: float, charge: int) -> float:
    ef = EmpiricalFormula(formula)
    mono_weight = ef.getMonoWeight()
    mz = (mono_weight + adduct_mass) / abs(charge)
    return mz

def get_element_counts(formula):
    counts = {}
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    for element, count in elements:
        counts[element] = int(count) if count else 1
    return counts

def simulate_variations(element_counts, variations):
    all_changes = {}
    for element, (min_change, max_change) in variations.items():
        all_changes[element] = range(min_change, max_change + 1)

    simulated_formulas = set()  # Use a set to avoid duplicates
    for changes in product(*all_changes.values()):
        new_counts = element_counts.copy()
        for i, element in enumerate(all_changes):
            new_counts[element] = max(0, new_counts.get(element, 0) + changes[i])
        new_formula = ''.join(f"{elem}{count}" for elem, count in new_counts.items() if count)
        simulated_formulas.add(new_formula)  # Add to set
    
    return simulated_formulas

def calculate_mz_from_file(filename, output_filename, adduct, mode, variations):
    results = []

    adduct_dict = adducts.positive_adducts if mode == 'pos' else adducts.negative_adducts

    with open(filename, 'r') as f:
        lines = f.readlines()

    for line in lines:
        smiles, substrate_name = line.strip().split()
        mol = Chem.MolFromSmiles(smiles)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        element_counts = get_element_counts(formula)
        simulated_formulas = simulate_variations(element_counts, variations)
        
        for sim_formula in simulated_formulas:
            adduct_mass, charge = adduct_dict[adduct]
            mz = calculateMass(sim_formula, adduct_mass, charge)
            results.append((substrate_name, sim_formula, mz))

    with open(output_filename, 'w') as f:
        f.write("Substrate,SimulatedFormula,Simulated-m/z\n")
        for substrate, formula, mz in results:
            f.write(f"{substrate},{formula},{mz}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate m/z of substrates from SMILES data.')
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input .smi file')
    parser.add_argument('--output_file', type=str, required=True, help='Path to save the output results')
    parser.add_argument('--parameters_file', type=str, required=True, help='Path to the variations.param file')

    args = parser.parse_args()

    # Load the parameters from the variations.param file
    with open(args.parameters_file, 'r') as file:
        params = yaml.safe_load(file)

    variations = {
        'H': params['H_variation'],
        'C': params['C_variation'],
        'O': params['O_variation'],
        'N': params['N_variation'],
        'F': params['F_variation'],
        'Cl': params['Cl_variation'],
    }

    calculate_mz_from_file(
        args.input_file, 
        args.output_file, 
        params['adduct'], 
        params['mode'], 
        variations
    )
