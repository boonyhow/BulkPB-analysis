#!/bin/bash

# Exit on error
set -e

ENV_NAME=$(grep '^name:' env_config_files/environment.yml | cut -d' ' -f2)
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

echo ">>> Running data download and preprocessing scripts..."
python3 code_modules/1_1_1_data_download.py
python3 code_modules/1_1_2_gencode_download.py
python3 code_modules/1_2_1_GSE176078_preprocessing.py
python3 code_modules/1_2_2_GSE217517_preprocessing.py
python3 code_modules/1_2_3_GSE226163_preprocessing.py

echo ">>> All scripts completed successfully."
