# Start from a clean Ubuntu image
FROM ubuntu:22.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies required by pipeline and Conda
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists\*

# Copy project files into container's filesystem
# Create a directory called /app to hold the project.
WORKDIR /app
COPY . .

# Download and install Mambaforge
RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh" -O mambaforge.sh && \
    bash mambaforge.sh -b -p /opt/conda && \
    rm mambaforge.sh
ENV PATH="/opt/conda/bin:${PATH}"

# Create the conda environment from environment.yml
RUN mamba env create -f environment.yml

# Activate the conda environment for all subsequent commands
SHELL ["mamba", "run", "-n", "nanopore_analysis", "/bin/bash", "-c"]

# Environment is now built and ready.
ENTRYPOINT ["python", "src/run_pipeline.py"]