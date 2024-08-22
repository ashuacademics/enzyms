#!/usr/bin/env python3
import argparse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from pyopenms import EmpiricalFormula
import re
from itertools import product

# Define dictionaries for positive and negative mode adducts with their respective mass shifts and charges
# The adducts are defined as a tuple of (mass_shift, charge)
# The mass shift is the difference between the mass of the adduct and the mass of the adduct ion
# Adduct masses are taken from https://fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator/

positive_adducts = {
    
    "[M+3H]3+": (1.007276, 3),
    "[M+2H+Na]3+": (8.334590, 3),
    "[M+H+2Na]3+": (15.7661904, 3),
    "[M+3Na]3+": (22.989218, 3),
    "[M+2H]2+": (1.007276, 2),
    "[M+H+NH4]2+": (9.520550, 2),
    "[M+H+Na]2+": (11.998247, 2),
	"[M+H+K]2+": (19.985217, 2),
	"[M+ACN+2H]2+": (21.520550, 2),
	"[M+2Na]2+": (22.989218, 2),
	"[M+2ACN+2H]2+": (42.033823, 2),
	"[M+3ACN+2H]2+": (62.547097, 2),
	"[M+H+]": (1.007276, 1),
	"[M+NH4]+": (18.033823, 1),
	"[M+Na]+": (22.989218, 1),
	"[M+CH3OH+H]+": (33.033489, 1),
	"[M+K]+": (38.963158, 1),
	"[M+ACN+H]+": (42.033823, 1),
	"[M+2Na-H]+": (44.971160, 1),
	"[M+IsoProp+H]+": (61.06534, 1),
	"[M+ACN+Na]+": (64.015765, 1),
	"[M+2K-H]+": (76.919040, 1),
	"[M+DMSO+H]+": (79.02122, 1),
	"[M+2ACN+H]+": (83.060370, 1),
	"[M+IsoProp+Na+H]+": (84.05511, 1),
	"[2M+H]+": (1.007276, 1),
	"[2M+NH4]+": (18.033823, 1),
	"[2M+Na]+": (22.989218, 1),
	"[2M+K]+": (38.963158, 1),
	"[2M+ACN+H]+": (42.033823, 1),
	"[2M+ACN+Na]+": (64.015765, 1),
       
}

negative_adducts = {
    
    "[M-3H]3-": (-1.007276, 3),
	"[M-2H]2-": (-1.007276, 2),
	"[M-H2O-H]-": (-19.01839, 1),
	"[M-H]-": (-1.007276, 1),
	"[M+Na-2H]-": (20.974666, 1),
	"[M+Cl]-": (34.969402, 1),
	"[M+K-2H]-": (36.948606, 1),
	"[M+FA-H]-": (44.998201, 1),
	"[M+Hac-H]-": (59.013851, 1),
	"[M+Br]-": (78.918885, 1),
	"[M+TFA-H]-": (112.985586, 1),
	"[2M-H]-": (-1.007276, 1),
	"[2M+FA-H]-": (44.998201, 1),
	"[2M+Hac-H]-": (59.013851, 1),
	"[3M-H]-": (-1.007276, 1),
    
	}

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

    adduct_dict = positive_adducts if mode == 'pos' else negative_adducts

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
    parser.add_argument('--adduct', type=str, required=True, help='The adduct ion type')
    parser.add_argument('--mode', type=str, choices=['pos', 'neg'], required=True, help='The ionization mode (positive or negative)')
    parser.add_argument('--H_variation', type=int, nargs=2, metavar=('MIN', 'MAX'), default=(-4, 4), help='Variation range for Hydrogen')
    parser.add_argument('--C_variation', type=int, nargs=2, metavar=('MIN', 'MAX'), default=(-1, 0), help='Variation range for Carbon')
    parser.add_argument('--O_variation', type=int, nargs=2, metavar=('MIN', 'MAX'), default=(-1, 3), help='Variation range for Oxygen')
    parser.add_argument('--N_variation', type=int, nargs=2, metavar=('MIN', 'MAX'), default=(0, 0), help='Variation range for Nitrogen')
    parser.add_argument('--F_variation', type=int, nargs=2, metavar=('MIN', 'MAX'), default=(-1, 0), help='Variation range for Fluorine')
    parser.add_argument('--Cl_variation', type=int, nargs=2, metavar=('MIN', 'MAX'), default=(0, 0), help='Variation range for Chlorine')

    args = parser.parse_args()

    variations = {
        'H': args.H_variation,
        'C': args.C_variation,
        'O': args.O_variation,
        'N': args.N_variation,
        'F': args.F_variation,
        'Cl': args.Cl_variation,
        
    }

    calculate_mz_from_file(args.input_file, args.output_file, args.adduct, args.mode, variations)
