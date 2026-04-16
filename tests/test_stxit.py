#!/usr/bin/env python3
"""
STXit Unit Tests
Placeholder test file for CI validation
"""

import unittest
import tempfile
import os
from pathlib import Path

# Import STXit modules
from stxit.stx_db import STXDatabase
from stxit.io_utils import validate_fasta_file, parse_query_genomes

class TestSTXDatabase(unittest.TestCase):
    """Test STX database functionality"""
    
    def test_database_loading(self):
        """Test STX database loads without errors"""
        try:
            db = STXDatabase()
            self.assertIsNotNone(db.sequences)
            self.assertIsNotNone(db.metadata)
        except Exception as e:
            self.fail(f"STX database failed to load: {e}")
    
    def test_database_validation(self):
        """Test STX database validation"""
        db = STXDatabase()
        issues = db.validate_database()
        
        # Should have minimal issues with bundled database
        critical_issues = [issue for issue in issues if 'Missing' in issue]
        self.assertEqual(len(critical_issues), 0, f"Critical database issues: {critical_issues}")
    
    def test_subtype_coverage(self):
        """Test that all expected subtypes are present"""
        db = STXDatabase()
        subtypes = db.list_subtypes()
        
        expected_subtypes = {
            'stx1a', 'stx1c', 'stx1d',
            'stx2a', 'stx2b', 'stx2c', 'stx2d', 'stx2e', 'stx2f', 'stx2g', 'stx2h'
        }
        
        self.assertTrue(expected_subtypes.issubset(subtypes), 
                       f"Missing subtypes: {expected_subtypes - subtypes}")

class TestIOUtils(unittest.TestCase):
    """Test I/O utility functions"""
    
    def test_fasta_validation(self):
        """Test FASTA file validation"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">test_seq\nATCGATCGATCG\n")
            temp_file = f.name
        
        try:
            is_valid, message = validate_fasta_file(temp_file)
            self.assertTrue(is_valid, f"Valid FASTA rejected: {message}")
        finally:
            os.unlink(temp_file)
    
    def test_query_parsing(self):
        """Test query genome string parsing"""
        # Single file
        result = parse_query_genomes("test.fasta")
        self.assertEqual(result, [("test.fasta", "test")])
        
        # File with strain name
        result = parse_query_genomes("test.fasta:EC4115")
        self.assertEqual(result, [("test.fasta", "EC4115")])
        
        # Multiple files
        result = parse_query_genomes("file1.fasta:strain1,file2.fasta:strain2")
        expected = [("file1.fasta", "strain1"), ("file2.fasta", "strain2")]
        self.assertEqual(result, expected)

class TestCLIIntegration(unittest.TestCase):
    """Test CLI integration"""
    
    def test_cli_import(self):
        """Test that CLI module imports without errors"""
        try:
            from stxit.cli import create_parser
            parser = create_parser()
            self.assertIsNotNone(parser)
        except Exception as e:
            self.fail(f"CLI import failed: {e}")
    
    def test_help_generation(self):
        """Test that help can be generated"""
        from stxit.cli import create_parser
        parser = create_parser()
        
        # Should not raise exception
        help_text = parser.format_help()
        self.assertIn('STXit', help_text)
        self.assertIn('--genome', help_text)

if __name__ == '__main__':
    unittest.main()
