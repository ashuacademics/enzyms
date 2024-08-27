from flask import Flask, request, render_template, redirect, url_for, send_file, send_from_directory, flash
import os
import subprocess
import pandas as pd
import zipfile
import csv

app = Flask(__name__)
app.config['SECRET_KEY'] = 'your_secret_key_here'
app.config['UPLOAD_FOLDER'] = 'uploads/'
app.config['OUTPUT_FOLDER'] = 'output/'
app.config['ZIP_FOLDER'] = 'zips/'
app.config['MZML_FOLDER'] = 'uploads/mzML-files/'

os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['OUTPUT_FOLDER'], exist_ok=True)
os.makedirs(app.config['ZIP_FOLDER'], exist_ok=True)
os.makedirs(app.config['MZML_FOLDER'], exist_ok=True)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/info')
def info():
    return render_template('info.html')

@app.route('/upload', methods=['POST'])
def upload_files():
    # Save uploaded files
    smiles_file = request.files['smiles_file']
    variations_file = request.files['variations_file']
    params_file = request.files['params_file']
    samples_file = request.files['samples_file']
    mzml_files = request.files.getlist('mzml_files')
    
    smiles_file_path = os.path.join(app.config['UPLOAD_FOLDER'], smiles_file.filename)
    variations_file_path = os.path.join(app.config['UPLOAD_FOLDER'], variations_file.filename)
    params_file_path = os.path.join(app.config['UPLOAD_FOLDER'], params_file.filename)
    samples_file_path = os.path.join(app.config['UPLOAD_FOLDER'], samples_file.filename)
    
    smiles_file.save(smiles_file_path)
    variations_file.save(variations_file_path)
    params_file.save(params_file_path)
    samples_file.save(samples_file_path)

    # Save mzML files
    for mzml_file in mzml_files:
        mzml_file_path = os.path.join(app.config['MZML_FOLDER'], mzml_file.filename)
        mzml_file.save(mzml_file_path)

    # Run Docker command
    docker_command = [
        "docker", "run", "--rm",
        "-v", f"{os.path.abspath(app.config['UPLOAD_FOLDER'])}:/usr/src/app/input",
        "-v", f"{os.path.abspath(app.config['OUTPUT_FOLDER'])}:/usr/src/app/output",
        "-v", f"{os.path.abspath(app.config['MZML_FOLDER'])}:/usr/src/app/input/mzML-files",
        "enzyms",
        "--smi_file", f"/usr/src/app/input/{smiles_file.filename}",
        "--variations_file", f"/usr/src/app/input/{variations_file.filename}",
        "--params_file", f"/usr/src/app/input/{params_file.filename}",
        "--samples_file", f"/usr/src/app/input/{samples_file.filename}",
        "--mzml_dir", "/usr/src/app/input/mzML-files"
    ]
    
    # Change to uploads directory and run the Docker container
    subprocess.run(docker_command, cwd=app.config['UPLOAD_FOLDER'])

    # Create a ZIP file of the output
    zip_filename = "output_files.zip"
    zip_filepath = os.path.join(app.config['ZIP_FOLDER'], zip_filename)
    with zipfile.ZipFile(zip_filepath, 'w') as zipf:
        for root, dirs, files in os.walk(app.config['OUTPUT_FOLDER']):
            for file in files:
                zipf.write(os.path.join(root, file), arcname=file)

    return redirect(url_for('results'))

def detect_delimiter(csv_file):
    with open(csv_file, 'r') as f:
        first_line = f.readline()
        sniffer = csv.Sniffer()
        delimiter = sniffer.sniff(first_line).delimiter
        return delimiter

@app.route('/results', methods=['GET'])
def results():
    # Display available files in the output directory
    files = os.listdir(app.config['OUTPUT_FOLDER'])
    return render_template('results.html', files=files)

@app.route('/output/<filename>')
def output_file(filename):
    # Serve images and CSV files
    if filename.endswith('.png') or filename.endswith('.jpg'):
        return send_from_directory(app.config['OUTPUT_FOLDER'], filename)
    elif filename.endswith('.csv'):
        return redirect(url_for('display_csv', filename=filename))
    else:
        return send_from_directory(app.config['OUTPUT_FOLDER'], filename)

@app.route('/csv/<filename>')
def display_csv(filename):
    csv_path = os.path.join(app.config['OUTPUT_FOLDER'], filename)
    delimiter = detect_delimiter(csv_path)
    df = pd.read_csv(csv_path, delimiter=delimiter)

    # Clean the data
    df.columns = df.columns.str.strip()
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    df.dropna(how='all', inplace=True)
    df = df[~(df.applymap(lambda x: isinstance(x, str) and x.strip() == '').all(axis=1))]
    pd.options.display.float_format = '{:.2f}'.format
    df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    # Generate the table HTML
    table_html = df.to_html(classes='data table table-striped', header="true", index=False)

    # Ensure no extraneous whitespace or characters
    table_html = table_html.strip()

    return render_template('display_csv.html', tables=table_html, filename=filename.rsplit('.', 1)[0])

@app.route('/download/<filename>')
def download_zip(filename):
    zip_path = os.path.join(app.config['ZIP_FOLDER'], filename)
    if os.path.exists(zip_path):
        return send_file(zip_path, as_attachment=True)
    else:
        flash('File not found!', 'danger')
        return redirect(url_for('results'))

if __name__ == '__main__':
    app.run(debug=True)
