#!/usr/bin/env python3
"""
STXit Setup Script
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read version from __init__.py
def get_version():
    init_file = Path(__file__).parent / "stxit" / "__init__.py"
    for line in init_file.read_text().split('\n'):
        if line.startswith('__version__'):
            return line.split('"')[1]
    return "1.0.3"

# Read long description from README
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text() if readme_file.exists() else ""

setup(
    name="stxit",
    version=get_version(),
    description="STXit - Shiga Toxin Detection and Analysis Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Bioinformatics Team",
    author_email="stxit@tool.dev",
    url="https://github.com/kramppe/STXit",
    license="MIT",
    
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'stxit': ['../databases/*']
    },
    
    python_requires=">=3.8",
    
    install_requires=[
        "biopython>=1.79",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "requests>=2.25.0",
        "pyyaml>=6.0",
        "click>=8.0.0",
    ],
    
    extras_require={
        'dev': [
            'pytest>=7.0.0',
            'pytest-cov>=3.0.0',
            'black>=22.0.0',
            'flake8>=4.0.0',
            'mypy>=0.900',
        ],
        'plotting': [
            'plotly>=5.0.0',
            'kaleido>=0.2.0',
        ]
    },
    
    entry_points={
        'console_scripts': [
            'stxit=stxit.cli:main',
        ],
    },
    
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],
    
    keywords="bioinformatics genomics shiga-toxin stx prophage blast",
    
    project_urls={
        "Bug Reports": "https://github.com/kramppe/STXit/issues",
        "Source": "https://github.com/kramppe/STXit",
        "Documentation": "https://github.com/kramppe/STXit/blob/main/README.md",
    },
    
    zip_safe=False,
)
