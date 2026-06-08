# Use a pinned Python runtime to keep Docker rebuilds reproducible.
FROM python:3.11.15

# Set the working directory inside the container
WORKDIR /usr/src/app

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Install dependencies
RUN pip install --no-cache-dir \
    asari-metabolomics==1.17.1 \
    rdkit==2026.3.3 \
    pyopenms==3.2.0 \
    pandas==2.3.3 \
    matplotlib==3.10.9 \
    numpy==1.24.4 \
    PyYAML==6.0.3

# Make the script executable
RUN chmod +x run_pipeline.sh

# Set the default command to run the script with the provided arguments
ENTRYPOINT ["./run_pipeline.sh"]
