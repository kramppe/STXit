#!/bin/bash
set -euo pipefail

# STXit Installation Script
# Installs conda environment and STXit package

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "STXit Installation Script"
echo "========================="

# Check if conda/mamba is available
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "Using mamba for package management"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    echo "Using conda for package management"
else
    echo "Error: conda or mamba not found. Please install Miniconda or Mambaforge first."
    echo "Download from: https://github.com/conda-forge/miniforge#mambaforge"
    exit 1
fi

# Create conda environment
echo ""
echo "Creating STXit conda environment..."
if $CONDA_CMD env list | grep -q "^stxit "; then
    echo "Environment 'stxit' already exists."
    read -p "Remove and recreate? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        $CONDA_CMD env remove -n stxit -y
        $CONDA_CMD env create -f environment.yml
    else
        echo "Using existing environment."
    fi
else
    $CONDA_CMD env create -f environment.yml
fi

# Activate environment and install STXit
echo ""
echo "Installing STXit package..."

# Source conda activation
CONDA_BASE=$($CONDA_CMD info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate environment
conda activate stxit

# Install STXit in development mode
pip install -e .

echo ""
echo "Validating installation..."

# Test basic functionality
stxit --version
stxit --help > /dev/null

# Validate STX database
python -c "
from stxit.stx_db import STXDatabase
db = STXDatabase()
issues = db.validate_database()
if issues:
    print('Database validation issues:', issues)
    exit(1)
print('✓ STX database validation passed')
"

echo ""
echo "✓ STXit installation completed successfully!"
echo ""
echo "Usage:"
echo "  conda activate stxit"
echo "  stxit --genome genome.fasta --output results/"
echo ""
echo "Test with:"
echo "  # Download test genome"
echo "  wget -O EC4115.fasta 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_011353.1&rettype=fasta'"
echo "  # Run analysis"
echo "  stxit --genome EC4115.fasta --output test_results/"
echo ""
echo "Documentation: https://github.com/kramppe/STXit"
