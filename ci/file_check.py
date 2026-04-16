#!/usr/bin/env python3
"""
CI File Checker
Validates presence of required project files
"""

import os
import sys
from pathlib import Path

REQUIRED_FILES = [
    'README.md',
    'environment.yml',
    'setup.py',
    'Dockerfile',
    '.devcontainer/devcontainer.json',
    'stxit/__init__.py',
    'stxit/cli.py',
    'stxit/stx_typing.py',
    'stxit/stx_db.py',
    'stxit/reporting.py',
    'stxit/plotting.py',
    'stxit/variant_analysis.py',
    'stxit/phastest_client.py',
    'stxit/trna_runner.py',
    'stxit/io_utils.py',
    'stxit/phage_prediction.py',
    'databases/stx_references.fasta',
    'databases/stx_references.tsv',
    'pipeline_diagram.svg',
    '.github/workflows/ci.yml'
]

OPTIONAL_FILES = [
    'LICENSE',
    'CHANGELOG.md',
    'docs/',
    'tests/',
    'examples/'
]

def check_required_files():
    """Check that all required files exist"""
    
    print(f"Checking {len(REQUIRED_FILES)} required files...")
    
    missing = []
    present = []
    
    for file_path in REQUIRED_FILES:
        if os.path.exists(file_path):
            present.append(file_path)
            print(f"✓ {file_path}")
        else:
            missing.append(file_path)
            print(f"✗ {file_path} (missing)")
    
    print(f"\nOptional files:")
    for file_path in OPTIONAL_FILES:
        if os.path.exists(file_path):
            print(f"✓ {file_path}")
        else:
            print(f"- {file_path} (optional)")
    
    if missing:
        print(f"\n✗ {len(missing)} required files missing:")
        for file_path in missing:
            print(f"  - {file_path}")
        return False
    else:
        print(f"\n✓ All {len(REQUIRED_FILES)} required files present")
        return True

if __name__ == '__main__':
    success = check_required_files()
    sys.exit(0 if success else 1)
