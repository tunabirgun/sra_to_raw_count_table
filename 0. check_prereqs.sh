#!/bin/bash

# List of required Conda packages
REQUIRED_PACKAGES=("fastqc" "fastp" "multiqc" "star" "samtools" "subread")

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    echo -e "\033[0;31mError: Conda is not installed. Please install Miniconda or Anaconda first.\033[0m"
    exit 1
fi

# Environment check
echo -e "\033[0;34mWould you like to create a new Conda environment for RNA-seq analysis? (yes/no)\033[0m"
read -r create_env

if [[ "$create_env" == "yes" ]]; then
    echo -e "\033[0;34mEnter the name of the new Conda environment:\033[0m"
    read -r ENV_NAME

    # Create and activate the new environment
    echo -e "\033[0;34mCreating Conda environment: $ENV_NAME\033[0m"
    conda create -y -n "$ENV_NAME" fastqc fastp multiqc star samtools subread

    echo -e "\033[0;34mActivating Conda environment: $ENV_NAME\033[0m"
    source activate "$ENV_NAME"

elif [[ "$create_env" == "no" ]]; then
    echo -e "\033[0;34mProceeding with installation in the current Conda environment.\033[0m"
else
    echo -e "\033[0;31mInvalid response. Please run the script again and enter 'yes' or 'no'.\033[0m"
    exit 1
fi

# Function to check and install missing predefined tools/packages
check_install() {
    local package="$1"
    if ! conda list | grep -q "^$package"; then
        echo -e "\033[0;33mInstalling missing package: $package\033[0m"
        conda install -y -c bioconda "$package"
    else
        echo -e "\033[0;32m$package is already installed.\033[0m"
    fi
}

# Check and install each package
for package in "${REQUIRED_PACKAGES[@]}"; do
    check_install "$package"
done

echo -e "\033[0;32mAll dependencies are installed and up to date!\033[0m"
