FROM gbuck10/nanopore-analysis-base:latest

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Copy project files into container's filesystem
# Create a directory called /app to hold the project.
WORKDIR /app
ENV PYTHONPATH="/app"

# Add metadata to the image
LABEL version="1.0.0"

# Activate the conda environment for all subsequent commands
SHELL ["mamba", "run", "-n", "nanopore_analysis", "/bin/bash", "-c"]
COPY . .

RUN chmod +x externals/wgbs_tools/wgbstools externals/UXM_deconv/uxm && \
    ln -s /app/externals/wgbs_tools/wgbstools /opt/conda/envs/nanopore_analysis/bin/wgbstools && \
    ln -s /app/externals/UXM_deconv/uxm /opt/conda/envs/nanopore_analysis/bin/uxm

# Install the pipeline entry point
RUN pip install -e .


# Download and install Dorado
RUN mkdir -p tools && \
    version=$(yq '.pipeline_steps.setup.params.dorado_version' config.yaml) && \
    wget "https://cdn.oxfordnanoportal.com/software/analysis/dorado-${version}-linux-x64.tar.gz" -O tools/dorado.tar.gz && \
    tar -xzf tools/dorado.tar.gz -C tools/ && \
    rm tools/dorado.tar.gz && \
    echo "Dorado ${version} installed successfully"

# Write runtime_config.yaml with Dorado path
RUN version=$(yq '.pipeline_steps.setup.params.dorado_version' config.yaml) && \
    echo "tools:" > runtime_config.yaml && \
    echo "  dorado: /app/tools/dorado-${version}-linux-x64/bin/dorado" >> runtime_config.yaml

# Verify the installation
RUN nanopore_analysis --help

# Environment is now built and ready.
ENTRYPOINT ["nanopore_analysis"]
CMD ["--help"]