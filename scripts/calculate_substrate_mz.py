import argparse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from pyopenms import EmpiricalFormula

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
    mono_mass = ef.getMonoWeight()
    mz = (mono_mass + adduct_mass) / abs(charge)
    return mz

def calculate_mz_from_file(filename, output_filename, adduct, mode):
    results = []

    # Select the correct adducts dictionary based on the mode
    adduct_dict = positive_adducts if mode == 'pos' else negative_adducts

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
        formatter_class=argparse.RawDescriptionHelpFormatter  # This will respect the line breaks in the epilog string
    )
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input .smi file')
    parser.add_argument('--output_file', type=str, required=True, help='Path to save the output results')
    parser.add_argument('--adduct', type=str, required=True, help='The adduct ion type in quotes (e.g. "[M+H+]")')
    parser.add_argument('--mode', type=str, choices=['pos', 'neg'], required=True, help='The ionization mode (pos for positive, neg for negative)')
    
    args = parser.parse_args()
    calculate_mz_from_file(args.input_file, args.output_file, args.adduct, args.mode)
