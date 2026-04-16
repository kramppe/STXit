#!/usr/bin/env python3
"""
CI YAML Checker
Validates YAML configuration files
"""

import yaml
import sys
import os

def check_yaml_files():
    """Check YAML files for syntax errors"""
    
    yaml_files = [
        'environment.yml',
        '.github/workflows/ci.yml',
        '.devcontainer/devcontainer.json'  # JSON but we can check syntax
    ]
    
    print(f"Checking {len(yaml_files)} configuration files...")
    
    errors = []
    
    for file_path in yaml_files:
        if not os.path.exists(file_path):
            error_msg = f"Configuration file missing: {file_path}"
            print(f"✗ {error_msg}")
            errors.append(error_msg)
            continue
        
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            if file_path.endswith('.json'):
                import json
                json.loads(content)
            else:
                yaml.safe_load(content)
            
            print(f"✓ {file_path}")
            
        except yaml.YAMLError as e:
            error_msg = f"YAML syntax error in {file_path}: {e}"
            print(f"✗ {error_msg}")
            errors.append(error_msg)
            
        except json.JSONDecodeError as e:
            error_msg = f"JSON syntax error in {file_path}: {e}"
            print(f"✗ {error_msg}")
            errors.append(error_msg)
            
        except Exception as e:
            error_msg = f"Error checking {file_path}: {e}"
            print(f"✗ {error_msg}")
            errors.append(error_msg)
    
    # Validate environment.yml structure
    if os.path.exists('environment.yml'):
        try:
            with open('environment.yml', 'r') as f:
                env_data = yaml.safe_load(f)
            
            required_keys = ['name', 'channels', 'dependencies']
            missing_keys = [key for key in required_keys if key not in env_data]
            
            if missing_keys:
                error_msg = f"environment.yml missing keys: {missing_keys}"
                print(f"✗ {error_msg}")
                errors.append(error_msg)
            else:
                print(f"✓ environment.yml structure valid")
                
                # Check for essential dependencies
                deps = env_data.get('dependencies', [])
                essential = ['python', 'blast', 'biopython']
                
                missing_deps = []
                for dep in essential:
                    if not any(str(d).startswith(dep) for d in deps):
                        missing_deps.append(dep)
                
                if missing_deps:
                    error_msg = f"environment.yml missing essential dependencies: {missing_deps}"
                    print(f"⚠ {error_msg}")
                    # Warning, not error
                else:
                    print(f"✓ environment.yml essential dependencies present")
        
        except Exception as e:
            error_msg = f"Error validating environment.yml structure: {e}"
            print(f"✗ {error_msg}")
            errors.append(error_msg)
    
    if errors:
        print(f"\n✗ {len(errors)} configuration errors found:")
        for error in errors:
            print(f"  - {error}")
        return False
    else:
        print(f"\n✓ All configuration files valid")
        return True

if __name__ == '__main__':
    success = check_yaml_files()
    sys.exit(0 if success else 1)
