# Use an official Python runtime as a parent image
FROM python:3.11

# Set the working directory inside the container
WORKDIR /usr/src/app

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Install dependencies
RUN pip install --no-cache-dir \
    asari-metabolomics \
    rdkit \
    pyopenms \
    pandas \
    matplotlib \
    numpy==1.24.4

# Make the script executable
RUN chmod +x run_pipeline.sh

# Set the default command to run the script with the provided arguments
ENTRYPOINT ["./run_pipeline.sh"]
