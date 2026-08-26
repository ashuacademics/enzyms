from flask import Flask, request, render_template, redirect, url_for, send_file, send_from_directory, flash
import os
import subprocess
import pandas as pd
import zipfile
import csv
import uuid

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
    # Generate a unique ID for this upload session
    session_id = str(uuid.uuid4())
    
    # Create session-specific directories
    session_upload_folder = os.path.join(app.config['UPLOAD_FOLDER'], session_id)
    session_output_folder = os.path.join(app.config['OUTPUT_FOLDER'], session_id)
    session_mzml_folder = os.path.join(session_upload_folder, 'mzML-files')

    os.makedirs(session_upload_folder, exist_ok=True)
    os.makedirs(session_output_folder, exist_ok=True)
    os.makedirs(session_mzml_folder, exist_ok=True)
        
    # Save uploaded files
    smiles_file = request.files['smiles_file']
    variations_file = request.files['variations_file']
    params_file = request.files['params_file']
    samples_file = request.files['samples_file']
    mzml_files = request.files.getlist('mzml_files')
    
    smiles_file.save(os.path.join(session_upload_folder, smiles_file.filename))
    variations_file.save(os.path.join(session_upload_folder, variations_file.filename))
    params_file.save(os.path.join(session_upload_folder, params_file.filename))
    samples_file.save(os.path.join(session_upload_folder, samples_file.filename))

    # Save mzML files
    for mzml_file in mzml_files:
        mzml_file.save(os.path.join(session_mzml_folder, mzml_file.filename))

    # Run Docker command
    docker_command = [
        "docker", "run", "--rm",
        "-v", f"{os.path.abspath(session_upload_folder)}:/usr/src/app/input",
        "-v", f"{os.path.abspath(session_output_folder)}:/usr/src/app/output",
        "-v", f"{os.path.abspath(session_mzml_folder)}:/usr/src/app/input/mzML-files",
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
    zip_filename = f"{session_id}.zip"
    zip_filepath = os.path.join(app.config['ZIP_FOLDER'], zip_filename)
    with zipfile.ZipFile(zip_filepath, 'w') as zipf:
        for root, dirs, files in os.walk(app.config['OUTPUT_FOLDER']):
            for file in files:
                zipf.write(os.path.join(root, file), arcname=file)

    return redirect(url_for('results', session_id=session_id))

def detect_delimiter(csv_file):
    with open(csv_file, 'r') as f:
        first_line = f.readline()
        sniffer = csv.Sniffer()
        delimiter = sniffer.sniff(first_line).delimiter
        return delimiter

@app.route('/results/<session_id>', methods=['GET'])
def results(session_id):
    # Display available files in the output directory
    session_output_folder = os.path.join(app.config['OUTPUT_FOLDER'], session_id)
    files = os.listdir(session_output_folder)
    return render_template('results.html', files=files, session_id=session_id)

@app.route('/output/<session_id>/<filename>')
def output_file(session_id, filename):
    # Serve images and CSV files
    session_output_folder = os.path.join(app.config['OUTPUT_FOLDER'], session_id)
    if filename.endswith('.csv'):
        return redirect(url_for('display_csv', session_id=session_id, filename=filename))
    return send_from_directory(session_output_folder, filename)

@app.route('/csv/<session_id>/<filename>')
def display_csv(session_id, filename):
    csv_path = os.path.join(app.config['OUTPUT_FOLDER'], session_id, filename)
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
    app.run(host="0.0.0.0", port=5000)