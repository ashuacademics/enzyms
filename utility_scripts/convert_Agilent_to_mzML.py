import os
import subprocess
import argparse

def convert_d_to_mzml(input_file, output_directory):
    """
    Use Dockerized msConvert (from chambm/pwiz-skyline-i-agree-to-the-vendor-licenses) 
    to convert .d to mzML format.
    """
    input_dir = os.path.dirname(os.path.abspath(input_file))
    input_filename = os.path.basename(input_file)
    
    cmd = [
        "docker", "run", "--rm",
        "-e", "WINEDEBUG=-all",
        "-v", f"{input_dir}:/data/input",
        "-v", f"{os.path.abspath(output_directory)}:/data/output",
        "chambm/pwiz-skyline-i-agree-to-the-vendor-licenses",
        "wine", "msconvert", f"/data/input/{input_filename}", "-o", "/data/output", "--mzML"
    ]
    subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser(description='Convert .d file to mzML format using Dockerized msConvert')
    parser.add_argument('--input_d_file', type=str, required=True, help='Path to the .d file')
    parser.add_argument('--output_directory', type=str, required=True, help='Directory where the mzML file will be saved')
    args = parser.parse_args()

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    convert_d_to_mzml(args.input_d_file, args.output_directory)
    print(f"Conversion to mzML finished. Output saved in: {args.output_directory}")

if __name__ == "__main__":
    main()
