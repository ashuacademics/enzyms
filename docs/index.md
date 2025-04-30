# EnzyMS: A LC-MS Data Analysis Pipeline for Enzyme Biocatalysis

**EnzyMS** is a modular and user-friendly pipeline designed to streamline the analysis of high-resolution LC-MS data from biocatalytic experiments. It is especially useful when testing the conversion of one or more substrates by an enzyme library (e.g., wildtype or engineered variants) in 96-well or high-throughput formats.

The workflow supports:

- Conversion of raw LC-MS data to open formats (mzML)  
- Detection and alignment of chromatographic features using [Asari](https://asari.readthedocs.io/)  
- Background subtraction using empty vector control (EVC) samples  
- Prediction and filtering of expected biotransformation products using [RDKit](https://www.rdkit.org/) and [pyOpenMS](https://pyopenms.readthedocs.io/) 
- Visualization of significant product peaks

The pipeline is modular and each step can be run independently. It is customizable for different enzyme classes and expected transformations.


