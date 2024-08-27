import yaml
import argparse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from pyopenms import EmpiricalFormula
import adducts  


'''
This script calculates the m/z of a substrate from its SMILES data. The script reads a .smi file containing the SMILES data of the substrates and their names. 
The script calculates the m/z of the substrates using the adduct ion type and the ionization mode provided in the variations.param file. The results are saved in a .csv file.

Usage: 

python calculate_substrate_mz.py --input_file substrates.smi --output_file substrate_mz.csv --parameters_file variations.param

Arguments:
--input_file: Path to the input .smi file
--output_file: Path to save the output results
--parameters_file: Path to the variations.param file

The variations.param file should include:
adduct: The adduct ion type in quotes (e.g. "[M+H+]")
mode: The ionization mode (pos for positive, neg for negative)
'''

def calculateMass(formula: str, adduct_mass: float, charge: int) -> float:
    ef = EmpiricalFormula(formula)
    mono_mass = ef.getMonoWeight()
    mz = (mono_mass + adduct_mass) / abs(charge)
    return mz

def calculate_mz_from_file(filename, output_filename, adduct, mode):
    results = []

    # Select the correct adducts dictionary based on the mode
    adduct_dict = adducts.positive_adducts if mode == 'pos' else adducts.negative_adducts

    with open(filename, 'r') as f:
        lines = f.readlines()

    for line in lines:
        smiles, substrate_name = line.strip().split()
        mol = Chem.MolFromSmiles(smiles)
        formula = rdMolDescriptors.CalcMolFormula(mol)

        adduct_mass, charge = adduct_dict[adduct]
        mz_substrate = calculateMass(formula, adduct_mass, charge)
        results.append((substrate_name, mz_substrate))

    with open(output_filename, 'w') as f:
        f.write("Substrate,Substrate-m/z\n")
        for res in results:
            f.write(f"{res[0]},{repr(res[1])}\n")


if __name__ == "__main__":
    # Examples of common adduct ions to be displayed after the help message
    adduct_examples = (
        "Examples of common adduct ions:\n"
        "Positive Mode: [M+H+], [M+Na]+, [M+K]+, [M+NH4]+, [M+2H]2+\n"
        "Negative Mode: [M-H]-, [M+Cl]-, [M-H2O-H]-, [2M-H]-, [M+FA-H]-"
    )

    parser = argparse.ArgumentParser(
        description='Calculate m/z of substrate from SMILES data.',
        epilog=adduct_examples,
        formatter_class=argparse.RawDescriptionHelpFormatter  
    )
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input .smi file')
    parser.add_argument('--output_file', type=str, required=True, help='Path to save the output results')
    parser.add_argument('--parameters_file', type=str, required=True, help='Path to the variations.param file')
    
    args = parser.parse_args()

    # Load the parameters from the variations.param file
    with open(args.parameters_file, 'r') as file:
        params = yaml.safe_load(file)

    calculate_mz_from_file(args.input_file, args.output_file, params['adduct'], params['mode'])
