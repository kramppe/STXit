#!/usr/bin/env python3
"""
CI Syntax Checker
Validates Python syntax for all .py files
"""

import ast
import sys
import os
from pathlib import Path

def check_python_syntax():
    """Check syntax of all Python files"""
    
    python_files = []
    
    # Find all Python files
    for root, dirs, files in os.walk('.'):
        # Skip hidden directories and __pycache__
        dirs[:] = [d for d in dirs if not d.startswith('.') and d != '__pycache__']
        
        for file in files:
            if file.endswith('.py'):
                python_files.append(os.path.join(root, file))
    
    print(f"Checking syntax of {len(python_files)} Python files...")
    
    errors = []
    
    for file_path in python_files:
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                source = f.read()
            
            # Parse syntax
            ast.parse(source, filename=file_path)
            print(f"✓ {file_path}")
            
        except SyntaxError as e:
            error_msg = f"Syntax error in {file_path}:{e.lineno}: {e.msg}"
            print(f"✗ {error_msg}")
            errors.append(error_msg)
            
        except UnicodeDecodeError as e:
            error_msg = f"Encoding error in {file_path}: {e}"
            print(f"✗ {error_msg}")
            errors.append(error_msg)
            
        except Exception as e:
            error_msg = f"Error checking {file_path}: {e}"
            print(f"✗ {error_msg}")
            errors.append(error_msg)
    
    if errors:
        print(f"\n{len(errors)} syntax errors found:")
        for error in errors:
            print(f"  - {error}")
        return False
    else:
        print(f"\n✓ All {len(python_files)} Python files passed syntax check")
        return True

if __name__ == '__main__':
    success = check_python_syntax()
    sys.exit(0 if success else 1)
