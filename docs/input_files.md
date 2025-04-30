# Required Input Files for EnzyMS

To run the EnzyMS pipeline successfully, you must provide the following five types of input files. Each plays a crucial role in detecting, matching, and visualizing molecular features from LCMS-QTOF data. All files must be placed in 'input' directory of EnzyMS root directory.

## 📄 1. SMILES File (`*.smi`)

- Filename example: `SoraphenA.smi`
- Contains the molecular structure of the target compound using [SMILES notation](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system).
- Format: One SMILES string with compound name per file.
- Example:
 
```
C[C@H]1/C=C/[C@H]([C@H](CCCC[C@H](OC(=O)[C@H]([C@@]2([C@@H]([C@H]([C@@H]([C@H]1O2)C)O)OC)O)C)C3=CC=CC=C3)OC)OC SoraphenA
```
## ⚙️ 2. Variations Parameter File (`variations.param`)

- Describes atomic-level variations and adducts to predict possible biotransformations.
- Required fields typically include:
    - Ionization mode (`pos` or `neg`)
    - Expected adduct ions (e.g., `[M+H]+`, `[M+Na]+`)
    - Atomic substitutions or modifications (e.g., `-CH3 +H`)
- This file is used to generate anticipated product m/z values.

- Filename example: `variations.param`
- Example:

```
adduct: "[M+Na]+"
mode: pos
H_variation: [-3, 3]
C_variation: [-1, 0]
O_variation: [-1, 3]
N_variation: [0, 0]
F_variation: [0, 0]
Cl_variation: [0, 1]
num_peaks: 100
```
## 🛠 3. ASARI Parameters File (`parameters.yaml`)

- Configuration file for the [Asari](https://asari.readthedocs.io/) feature detection engine.
- This file contains parameters for retention time windows, peak detection, feature alignment etc
- For details, please check [this article](https://www.nature.com/articles/s41467-023-39889-1).
- Should match the experimental setup for best results.

- Filename example: `parameters.yaml`
- Example:

```
project_name: asari_project
mz_tolerance_ppm: 10
mode: pos
# Below rarely modified
multicores: 4
outdir: output
reference: null
autoheight: true
database_mode: auto
gaussian_shape: 0.3
min_intensity_threshold: 1000
min_peak_height: 10000
min_timepoints: 6
rt_align_method: lowess
rt_align_on: true
rtime_tolerance: 50
signal_noise_ratio: 5
wlen: 25
```

## 🧪 4. Sample List File (`list_of_samples.txt`)

- Lists all LC-MS samples and their type (e.g., samples, evc).
- Filename example: `list_of_samples.txt`
- Format: CSV (with header)
- Example:
```
Enzyme,Sample
Soraphen_A1,sample
Soraphen_A2,sample
Soraphen_A3,sample
Soraphen_A4,evc
Soraphen_A5,sample
Soraphen_A6,sample
Soraphen_A7,sample
Soraphen_A8,evc
Soraphen_A9,sample
```

## 📈 5. LC-MS Data Files (`*.mzML`)

- LC-MS data in `.mzML` format.
- Each file must correspond to a sample listed in `list_of_samples.txt`.
- All files must be placed in 'mzML-files' directory
- For conversion of LC-MS files from vendor format e.g. Agilent `.d` to `.mzML` use [ProteoWizard's `msConvert`](https://proteowizard.sourceforge.io/).
- You can also use 'convert_Agilent_to_mzML.py' utility_script to multiple .d files. Run following commands in directory where `.d` files are located

```
mkdir mzML-files
for d in *.d ; do python3 convert_Agilent_to_mzML.py --input_d_file $d --output_directory ./mzML-files ; done
```