# LCMS Data Analysis Pipeline

The LCMS (Liquid Chromatography-Mass Spectrometry) Data Analysis Pipeline is designed to streamline and automate the analysis of LCMS data. This pipeline takes input in the form of SMILES files, parameter files, and sample lists, processes the data using specified configurations, and outputs analyzed data in a variety of formats including CSV and images.

## Features:

- **Automated LCMS Data Analysis**: Streamline the processing and analysis of LCMS data with a single pipeline.
- **Multiple Input Types**: Accepts SMILES files, parameter configuration files, and sample lists.
- **Customizable Analysis**: Adjust analysis parameters through the parameter file.
- **Versatile Output**: Generates a range of output files, including CSVs and images, which can be viewed or downloaded.

## Required Input

To run the LCMS Data Analysis Pipeline, you need the following inputs:

- **SMILES File** (`*.smi`): Contains the molecular structure in SMILES format.
- **Parameter File** (`parameters.yaml`): Configuration file specifying the analysis parameters.
- **Sample List File** (`samples.txt`): Contains the list of samples to be analyzed.
- **mzML Files** (`*.mzML`): Raw data files from LCMS to be analyzed.

## Installation and Usage

You can run the LCMS Data Analysis Pipeline in several ways:

### 1. Running as a Docker Container

#### Prerequisites
[Docker](https://docs.docker.com/get-docker/) must be installed on your system.

#### Steps to Run
1. **Build the Docker Image**:
    ```bash
    docker build -t lcms-pipeline .
    ```
2. **Run the Pipeline**:
    ```bash
    docker run --rm -v $(pwd)/input:/usr/src/app/input -v $(pwd)/output:/usr/src/app/output lcms-pipeline --smi_file /usr/src/app/input/CPD041.smi --params_file /usr/src/app/input/parameters.yaml --samples_file /usr/src/app/input/list_of_samples.txt --mzml_dir /usr/src/app/input/mzML-files
    ```
    - `input`: Directory containing the input files.
    - `output`: Directory where the results will be saved.

### 2. Running as a Conda Environment

#### Prerequisites
[Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) must be installed.

#### Steps to Install and Run
1. **Create a Conda Environment**:
    ```bash
    conda create -n lcms-pipeline python=3.8
    conda activate lcms-pipeline
    ```
2. **Install Required Packages**:
    ```bash
    pip install -r requirements.txt
    ```
3. **Run the Shell Script**:
    ```bash
    bash pipeline.sh --smi_file input/CPD041.smi --params_file input/parameters.yaml --samples_file input/list_of_samples.txt --mzml_dir input/mzML-files
    ```

### 3. Running as a Web Application

#### Prerequisites
- [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) must be installed.
- [Flask](https://flask.palletsprojects.com/) and other dependencies as listed in `requirements.txt`.

#### Steps to Run
1. **Create a Conda Environment**:
    ```bash
    conda create -n lcms-webapp python=3.8
    conda activate lcms-webapp
    ```
2. **Install Required Packages**:
    ```bash
    pip install -r requirements.txt
    ```
3. **Run the Flask Application**:
    ```bash
    python app.py
    ```
4. **Access the Web Application**:
    Open your web browser and navigate to [http://127.0.0.1:5000](http://127.0.0.1:5000). You can upload the required files and run the pipeline through the web interface.

## Example Run

### Running via Docker
Assuming your input files are located in the `input/` directory and you want to save the output to the `output/` directory:

```bash
docker run --rm -v $(pwd)/input:/usr/src/app/input -v $(pwd)/output:/usr/src/app/output lcms-pipeline --smi_file /usr/src/app/input/CPD041.smi --params_file /usr/src/app/input/parameters.yaml --samples_file /usr/src/app/input/list_of_samples.txt --mzml_dir /usr/src/app/input/mzML-files
