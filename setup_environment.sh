#!/bin/bash
# Setup script for hifi-human-wgs conda environment

set -e

echo "Setting up hifi-human-wgs conda environment..."

# Create conda environment from environment.yml
if conda env list | grep -q "hifi-human-wgs"; then
    echo "Environment 'hifi-human-wgs' already exists."
    read -p "Do you want to remove and recreate it? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n hifi-human-wgs -y
        echo "Creating new environment..."
        conda env create -f environment.yml
    else
        echo "Updating existing environment..."
        conda env update -n hifi-human-wgs -f environment.yml --prune
    fi
else
    echo "Creating new environment..."
    conda env create -f environment.yml
fi

echo ""
echo "Environment setup complete!"
echo ""
echo "To activate the environment, run:"
echo "    conda activate hifi-human-wgs"
echo ""
echo "To verify the installation, run:"
echo "    miniwdl --version"
echo ""
