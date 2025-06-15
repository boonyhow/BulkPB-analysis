#!/bin/bash

# Exit on error
set -e

echo ">>> Creating Conda environment from environment.yml..."
conda env create -f env_config_files/environment.yml -y

ENV_NAME=$(grep '^name:' env_config_files/environment.yml | cut -d' ' -f2)
echo ">>> Environment '$ENV_NAME' created."

echo ">>> Activating environment and installing R packages..."
# Activate env and run R script for package installation
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

Rscript code_modules/0_install_r_packages.R

echo ">>> Environment setup complete."
