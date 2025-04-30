# Installation & Running EnzyMS

The EnzyMS pipeline can be run in three ways:

- **Docker** (Recommended) — Easy, isolated, and portable
- **Conda Environment** — Ideal for local development and command-line usage
- **Local Web App (Flask GUI)** — Best for GUI-based workflows

---

## 🚀 Option 1: Run with Docker (Recommended)

### ✅ Prerequisites
- 📦 Install Docker: [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/)

### 🛠️ Clone EnzyMS from github repository
```
git clone https://github.com/ashuacademics/enzyms.git
cd enzyms
```

### 🔧 Build Docker Image 
```
docker build -t enzyms .
```

### 🧬 Run pipeline using docker command
```
docker run --rm \
  -v $(pwd)/input:/usr/src/app/input \                
  -v $(pwd)/output:/usr/src/app/output \              
  enzyms \                                            
  --smi_file /usr/src/app/input/SoraphenA.smi \        
  --variations_file /usr/src/app/input/variations.param \               
  --params_file /usr/src/app/input/parameters.yaml \   
  --samples_file /usr/src/app/input/list_of_samples.txt \ 
  --mzml_dir /usr/src/app/input/mzML-files             
```

## 🚀 Option 2: Run from a conda environment

### ✅  Prerequisites

- 📦 Install Anaconda: [https://www.anaconda.com/products/distribution](https://www.anaconda.com/products/distribution)

### 🛠️ Create a conda environment and install required packages

**Option 1: Create a fresh environment manually**
```
conda create -n enzyms python=3.11
conda activate enzyms
pip install asari-metabolomics rdkit pyopenms pandas matplotlib numpy==1.24.4
```

**Option 2: Create from a YAML environment file**
```
conda env create -f conda_environment.yaml
conda activate enzyms
```

### 🧬 Run pipeline using shell script
```
sh run_pipeline.sh \
    --smi_file input/SoraphenA.smi \
    --variations_file input/variations.param \
    --params_file input/parameters.yaml \
    --samples_file input/list_of_samples.txt \
    --mzml_dir input/mzML-files
```

## 🚀 Option 3: Run as a web application

### ✅ Prerequisites

- Install Anaconda: [https://www.anaconda.com/products/distribution](https://www.anaconda.com/products/distribution)

### 🛠️ Create a conda environment and install required packages

```
conda create -n enzyms-webapp python=3.11
conda activate enzyms-webapp
pip install Flask asari-metabolomics rdkit pyopenms pandas matplotlib numpy==1.24.4
```

### 🧬 Run pipeline's Flask application

```
python enzyms_app.py
```

### 🌐 Access the Web Application

- Open your web browser and go to: [http://127.0.0.1:5000](http://127.0.0.1:5000)  
- From the interface, you can upload the required files and run the EnzyMS pipeline directly through the web UI.
