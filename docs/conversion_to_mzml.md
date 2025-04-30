# Conversion of LC-MS Data Files (`*.mzML`)

To run the EnzyMS LC-MS Data Analysis Pipeline, input LC-MS data must be in the standard **`.mzML`** format.

- Each `.mzML` file should correspond to a sample listed in `list_of_samples.txt`.
- All converted files should be placed in the `mzML-files` directory.

---

## 📦 Conversion from Vendor Formats

To convert raw LC-MS files (e.g., Agilent `.d` directories) into `.mzML`, use [ProteoWizard's `msConvert`](https://proteowizard.sourceforge.io/). This tool supports a wide range of vendor formats, including:

- Agilent (`.d`)
- Thermo Scientific (`.RAW`)
- Shimadzu (`.LCD`)
- Waters (`.raw`)
- Bruker, AB Sciex, and many more

A full list of supported formats is available in the [official msConvert documentation](https://proteowizard.sourceforge.io/tools/msconvert.html).

> ⚠️ The EnzyMS pipeline has been tested only with **Agilent `.d`** files. However, it should work with any vendor data as long as it is successfully converted to `.mzML`.

---

## 🐳 Using Dockerized `msConvert` via a Python Script

We provide a utility script named `convert_Agilent_to_mzML.py`, which uses a [Docker container](https://hub.docker.com/repository/docker/proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses) for msConvert. This makes it easier to run on systems without native Windows installations of ProteoWizard.

### 🔁 Convert a Single `.d` File

```
python3 convert_Agilent_to_mzML.py \
    --input_d_file example.d \
    --output_directory ./mzML-files
```


### 🔁 Convert multiple `.d` Files

```
for d in *.d ; do
    python3 convert_Agilent_to_mzML.py \
        --input_d_file "$d" \
        --output_directory ./mzML-files
done
```