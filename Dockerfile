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

# Tells python to look for modules in the project root.
ENV PYTHONPATH="${PYTHONPATH}:/app"

# Download and install Miniforge
RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" -O miniforge.sh && \
    bash miniforge.sh -b -p /opt/conda && \
    rm miniforge.sh
ENV PATH="/opt/conda/bin:${PATH}"

# Copy environment.yml
COPY environment.yml .

# Create the conda environment from environment.yml
RUN mamba env create -f environment.yml

# Activate the conda environment for all subsequent commands
SHELL ["mamba", "run", "-n", "nanopore_analysis", "/bin/bash", "-c"]

COPY . .

RUN echo "Conda environment activated successfully" && \
    python --version && \
    yq --version

# Environment is now built and ready.
ENTRYPOINT ["mamba", "run", "-n", "nanopore_analysis", "python", "src/run_pipeline.py"]
CMD ["all"]