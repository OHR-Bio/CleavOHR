#!/usr/bin/env python3
"""
Protease Cleavage Site Analyzer

This script analyzes protein sequences for potential protease cleavage sites
based on a database of known cleavage patterns.

Author: Claude
Date: 2025
"""

import pandas as pd
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import re

class ProteaseAnalyzer:
    def __init__(self, database_file):
        """
        Initialize the analyzer with a protease cleavage database.
        
        Args:
            database_file: Path to the TSV file containing protease cleavage data
        """
        # Amino acid conversion dictionaries
        self.three_to_one = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
            'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
            'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
            'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
            # Handle some special cases that might appear in the database
            'TyI': 'Y',  # Iodinated tyrosine - treat as tyrosine
            'TyrI': 'Y',  # Alternative notation
            'Tyr(I2)': 'Y',  # Another notation for iodinated tyrosine
            'TyI2': 'Y',  # Another variant
            # Add more special cases as needed
        }
        
        self.one_to_three = {v: k for k, v in self.three_to_one.items()}
        
        self.database = self._load_database(database_file)
        self.cleavage_patterns = self._parse_cleavage_patterns()
    
    def _load_database(self, database_file, species="Escherichia coli"):
        """Load the protease cleavage database from TSV file."""
        # Try different encodings
        encodings_to_try = ['utf-8', 'latin-1', 'iso-8859-1', 'cp1252']
        
        for encoding in encodings_to_try:
            try:
                print(f"Trying to load database with {encoding} encoding...")
                
                # First, let's try with flexible column handling
                df = pd.read_csv(database_file, sep='\t', header=None, encoding=encoding,
                               on_bad_lines='skip', engine='python')
                
                print(f"Successfully loaded database with {encoding} encoding!")
                print(f"Database contains {len(df)} entries")

                break
                
            except UnicodeDecodeError:
                print(f"Failed with {encoding} encoding, trying next...")
                continue
            except Exception as e:
                print(f"Error with {encoding} encoding: {e}")
                continue
        else:
            # If all encodings fail, try reading with error handling
            try:
                print("Trying to load with error handling...")
                df = pd.read_csv(database_file, sep='\t', header=None, 
                               encoding='utf-8', encoding_errors='replace',
                               on_bad_lines='skip', engine='python')
                print(f"Loaded database with error handling - {len(df)} entries")
            except Exception as e:
                raise Exception(f"Error loading database with all methods: {e}")
        
        # Check the number of columns and adjust if necessary
        print(f"Database has {len(df.columns)} columns")
        
        # If we have more than 28 columns, we need to handle this
        if len(df.columns) >= 23:
            # Take only the first 23 columns that we need
            expected_columns = ['CLE_ID', 'Enzyme_ID', 'Substrate', 'Cleavage_Site', 
                              'P4', 'P3', 'P2', 'P1', 'P1_prime', 'P2_prime', 
                              'P3_prime', 'P4_prime', 'Reference', 'UniProt_ID', 
                              'Position', 'Organism', 'Protease_Name', 'Col_Q', 
                              'Col_R', 'Sequence_Range', 'Col_T', 'Col_U', 'Physiological']
            
            # Use only the columns we need
            df = df.iloc[:, :23]  # Take first 23 columns
            df.columns = expected_columns
        else:
            # Pad with None columns if we have fewer than expected
            while len(df.columns) < 23:
                df[f'Extra_Col_{len(df.columns)}'] = None
            
            expected_columns = ['CLE_ID', 'Enzyme_ID', 'Substrate', 'Cleavage_Site', 
                              'P4', 'P3', 'P2', 'P1', 'P1_prime', 'P2_prime', 
                              'P3_prime', 'P4_prime', 'Reference', 'UniProt_ID', 
                              'Position', 'Organism', 'Protease_Name', 'Col_Q', 
                              'Col_R', 'Sequence_Range', 'Col_T', 'Col_U', 'Physiological']
            
            df.columns = expected_columns[:len(df.columns)]
        
        # Clean up the data - remove quotes and handle NULL values
        for col in df.columns:
            if df[col].dtype == 'object':
                df[col] = df[col].astype(str).str.strip("'\"")
                df[col] = df[col].replace('NULL', None)
                df[col] = df[col].replace('nan', None)
                # Remove any replacement characters
                df[col] = df[col].str.replace('�', '', regex=False)
        
        # Print some statistics
        print(f"Final database shape: {df.shape}")

        # uniprot = pd.read_csv("db/uniprot.txt", header=None, encoding='latin-1', on_bad_lines='skip',encoding_errors='replace' )
        # uniprot.columns = ["id", "protein_id", "name", "species", "uniprot"]
        # print(f"Read uniprot successfully")

        # # keep only Escherichia coli results
        # uniprot = uniprot[uniprot['species'].str.contains(species, case=False, na=False)]
        # print(f"Filtered to {len(uniprot)} rows containing {species}")

        # df = df.merge(uniprot, left_on='Enzyme_ID', right_on='protein_id', how='inner')

        # print(f"Merged uniprot successfully")

        return df
    
    def _convert_to_single_letter(self, amino_acid):
        """
        Convert three-letter amino acid code to single letter code.
        
        Args:
            amino_acid: Three-letter amino acid code or special notation
            
        Returns:
            Single letter amino acid code or None if conversion fails
        """
        if pd.isna(amino_acid) or amino_acid == 'NULL' or amino_acid == '-':
            return amino_acid
            
        # Clean up the amino acid string
        aa = str(amino_acid).strip()
        
        # Handle special cases and direct mapping
        if aa in self.three_to_one:
            return self.three_to_one[aa]
        
        # Try to extract the base amino acid from modified forms
        # Handle cases like "Tyr(I2)" or "ac-Phe"
        if '(' in aa:
            base_aa = aa.split('(')[0]
            if base_aa in self.three_to_one:
                return self.three_to_one[base_aa]
        
        # Handle cases with prefixes like "ac-Phe"
        if '-' in aa:
            parts = aa.split('-')
            for part in parts:
                if part in self.three_to_one:
                    return self.three_to_one[part]
        
        # If it's already a single letter, return it
        if len(aa) == 1 and aa.isalpha():
            return aa.upper()
        
        # Log unknown amino acids for debugging
        #print(f"Warning: Unknown amino acid code: {aa}")
        return None

    def _parse_cleavage_patterns(self):
        """Parse cleavage patterns from the database."""
        patterns = []
        
        for _, row in self.database.iterrows():
            # Extract the scissile bond (P1 and P1') and convert to single letter
            p1_three = row['P1']
            p1_prime_three = row['P1_prime']
            
            p1 = self._convert_to_single_letter(p1_three)
            p1_prime = self._convert_to_single_letter(p1_prime_three)
            
            # Skip if either scissile bond position is missing or couldn't be converted
            if not p1 or not p1_prime or p1 == 'NULL' or p1_prime == 'NULL':
                continue
            
            # Extract neighboring amino acids and convert to single letter
            neighbors = {}
            neighbor_cols = ['P4', 'P3', 'P2', 'P2_prime', 'P3_prime', 'P4_prime']
            
            for col in neighbor_cols:
                three_letter = row[col]
                single_letter = self._convert_to_single_letter(three_letter)
                neighbors[col] = single_letter
            
            pattern = {
                'cle_id': row['CLE_ID'],
                'enzyme_id': row['Enzyme_ID'],
                'protease_name': row['Protease_Name'],
                'scissile_bond': (p1, p1_prime),
                'neighbors': neighbors,
                'substrate': row['Substrate'],
                'original_p1': p1_three,
                'original_p1_prime': p1_prime_three
            }
            
            patterns.append(pattern)
        
        return patterns
    
    def find_exact_scissile_matches(self, sequence):
        """
        Find exact matches for scissile bonds in the protein sequence.
        
        Args:
            sequence: Protein sequence string
            
        Returns:
            List of matches with protease information
        """
        matches = []
        sequence = str(sequence).upper()
        
        for pattern in self.cleavage_patterns:
            p1, p1_prime = pattern['scissile_bond']
            
            # Find all occurrences of the scissile bond
            for i in range(len(sequence) - 1):
                if sequence[i] == p1 and sequence[i + 1] == p1_prime:
                    match = {
                        'position': i + 1,  # 1-based position
                        'scissile_bond': f"{p1}|{p1_prime}",
                        'protease_name': pattern['protease_name'],
                        'enzyme_id': pattern['enzyme_id'],
                        'cle_id': pattern['cle_id'],
                        'match_type': 'exact_scissile',
                        'context': self._get_context(sequence, i, 3),
                        'original_p1': pattern.get('original_p1', p1),
                        'original_p1_prime': pattern.get('original_p1_prime', p1_prime)
                    }
                    matches.append(match)
        
        return matches
    
    def find_neighbor_matches(self, sequence, min_matches=3):
        """
        Find matches based on neighboring amino acids (at least 3/6 positions).
        
        Args:
            sequence: Protein sequence string
            min_matches: Minimum number of neighboring positions that must match
            
        Returns:
            List of matches with protease information
        """
        matches = []
        sequence = str(sequence).upper()
        
        for pattern in self.cleavage_patterns:
            p1, p1_prime = pattern['scissile_bond']
            
            # Find all occurrences of the scissile bond
            for i in range(len(sequence) - 1):
                if sequence[i] == p1 and sequence[i + 1] == p1_prime:
                    # Check neighboring amino acids
                    neighbor_matches = self._count_neighbor_matches(sequence, i, pattern['neighbors'])
                   
                    
                    if neighbor_matches >= min_matches:
                        match = {
                            'position': i + 1,  # 1-based position
                            'scissile_bond': f"{p1}|{p1_prime}",
                            'protease_name': pattern['protease_name'],
                            'enzyme_id': pattern['enzyme_id'],
                            'cle_id': pattern['cle_id'],
                            'match_type': 'neighbor_pattern',
                            'neighbor_matches': neighbor_matches,
                            'context': self._get_context(sequence, i, 5),
                            'original_p1': pattern.get('original_p1', p1),
                            'original_p1_prime': pattern.get('original_p1_prime', p1_prime)
                        }
                        matches.append(match)
        return matches
    
    def _count_neighbor_matches(self, sequence, scissile_pos, neighbors):
        """Count how many neighboring positions match the pattern."""
        matches = 0
        total_positions = 0
        
        # Define positions relative to scissile bond
        positions = {
            'P4': scissile_pos - 3,
            'P3': scissile_pos - 2,
            'P2': scissile_pos - 1,
            'P2_prime': scissile_pos + 2,
            'P3_prime': scissile_pos + 3,
            'P4_prime': scissile_pos + 4,
        }
        
        for pos_name, expected_aa in neighbors.items():
            # if pd.isna(expected_aa) or expected_aa == 'NULL' or expected_aa is None:
            #     continue
                
            seq_pos = positions[pos_name]
            
            # Check if position is within sequence bounds
            if 0 <= seq_pos < len(sequence):
                total_positions += 1
                
                # '-' means any amino acid is acceptable
                if expected_aa == '-' or sequence[seq_pos] == expected_aa: # or pd.isna(expected_aa) or expected_aa == 'NULL' or expected_aa is None:
                    matches += 1
        
        return matches
    
    def _get_context(self, sequence, scissile_pos, window=3):
        """Get the sequence context around the scissile bond."""
        start = max(0, scissile_pos - window)
        end = min(len(sequence), scissile_pos + window + 2)
        
        context = sequence[start:end]
        cleavage_site = scissile_pos - start
        
        # Insert cleavage indicator
        if cleavage_site < len(context) - 1:
            context = context[:cleavage_site + 1] + '|' + context[cleavage_site + 1:]
        
        return context
    
    def analyze_sequence(self, sequence, sequence_id="Unknown", min_neighbor_matches=3):
        """
        Analyze a protein sequence for cleavage sites.
        
        Args:
            sequence: Protein sequence string
            sequence_id: Identifier for the sequence
            min_neighbor_matches: Minimum neighbor matches required
            
        Returns:
            Dictionary with analysis results
        """
        # Find exact scissile matches
        exact_matches = self.find_exact_scissile_matches(sequence)
        
        # Find neighbor pattern matches
        neighbor_matches = self.find_neighbor_matches(sequence, min_neighbor_matches)
        
        # Remove duplicates (exact matches are also in neighbor matches)
        #unique_neighbor_matches = []
        #exact_positions = {(m['position'], m['cle_id']) for m in exact_matches}
        
        # for match in neighbor_matches:
        #     if (match['position'], match['cle_id']) not in exact_positions:
        #         unique_neighbor_matches.append(match)
        
        results = {
            'sequence_id': sequence_id,
            'sequence_length': len(sequence),
            'exact_matches': exact_matches,
            'neighbor_matches': neighbor_matches,
            'total_sites': len(exact_matches) + len(neighbor_matches)
        }
        
        return results
    
    def analyze_fasta(self, fasta_file, min_neighbor_matches=3):
        """
        Analyze all sequences in a FASTA file.
        
        Args:
            fasta_file: Path to FASTA file
            min_neighbor_matches: Minimum neighbor matches required
            
        Returns:
            List of analysis results for each sequence
        """
        results = []
        
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                result = self.analyze_sequence(
                    str(record.seq), 
                    record.id, 
                    min_neighbor_matches
                )
                results.append(result)
                
        except Exception as e:
            raise Exception(f"Error analyzing FASTA file: {e}")
        
        return results
    
    def print_results(self, results):
        """Print analysis results in a formatted way."""
        if isinstance(results, dict):
            results = [results]
        
        for result in results:
            print(f"\n{'='*60}")
            print(f"Sequence: {result['sequence_id']}")
            print(f"Length: {result['sequence_length']} amino acids")
            print(f"Total cleavage sites found: {result['total_sites']}")
            
            if result['exact_matches']:
                print(f"\nExact Scissile Bond Matches ({len(result['exact_matches'])}):")
                print("-" * 40)
                for match in result['exact_matches']:
                    print(f"Position {match['position']}: {match['scissile_bond']}")
                    print(f"  Protease: {match['protease_name']}")
                    print(f"  Enzyme ID: {match['enzyme_id']}")
                    print(f"  Context: {match['context']}")
                    print(f"  Database entry: {match.get('original_p1', '')}-{match.get('original_p1_prime', '')}")
                    print()
            
            if result['neighbor_matches']:
                print(f"\nNeighbor Pattern Matches ({len(result['neighbor_matches'])}):")
                print("-" * 40)
                for match in result['neighbor_matches']:
                    print(f"Position {match['position']}: {match['scissile_bond']}")
                    print(f"  Protease: {match['protease_name']}")
                    print(f"  Enzyme ID: {match['enzyme_id']}")
                    print(f"  Neighbor matches: {match['neighbor_matches']}/6")
                    print(f"  Context: {match['context']}")
                    print(f"  Database entry: {match.get('original_p1', '')}-{match.get('original_p1_prime', '')}")
                    print()
    
    def export_results(self, results, output_file):
        """Export results to a CSV file."""
        if isinstance(results, dict):
            results = [results]
        
        all_matches = []
        
        for result in results:
            # Add exact matches
            for match in result['exact_matches']:
                match_data = {
                    'sequence_id': result['sequence_id'],
                    'sequence_length': result['sequence_length'],
                    **match
                }
                all_matches.append(match_data)
            
            # Add neighbor matches
            for match in result['neighbor_matches']:
                match_data = {
                    'sequence_id': result['sequence_id'],
                    'sequence_length': result['sequence_length'],
                    **match
                }
                all_matches.append(match_data)
        
        df = pd.DataFrame(all_matches)
        df = df.drop_duplicates()  # Remove duplicate rows
        df.to_csv(output_file, index=False)
        print(f"Results exported to {output_file}")


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(description="Analyze protein sequences for protease cleavage sites")
    parser.add_argument("database", help="Path to the protease cleavage database (TSV file)")
    parser.add_argument("--sequence", "-s", help="Single protein sequence to analyze")
    parser.add_argument("--fasta", "-f", help="Path to FASTA file containing protein sequences")
    parser.add_argument("--min-neighbors", "-m", type=int, default=3, 
                       help="Minimum number of neighboring amino acids that must match (default: 3)")
    parser.add_argument("--output", "-o", help="Output CSV file for results")
    
    args = parser.parse_args()
    
    if not args.sequence and not args.fasta:
        print("Error: Please provide either a sequence (-s) or a FASTA file (-f)")
        return
    
    # Initialize analyzer
    analyzer = ProteaseAnalyzer(args.database)
    
    # Analyze sequences
    if args.sequence:
        results = analyzer.analyze_sequence(args.sequence, "Input_Sequence", args.min_neighbors)
    else:
        results = analyzer.analyze_fasta(args.fasta, args.min_neighbors)
    
    # Print results
    #analyzer.print_results(results)
    
    # Export results if requested
    if args.output:
        analyzer.export_results(results, args.output)


if __name__ == "__main__":
    main()


# Example usage as a module:
"""
# Initialize the analyzer
analyzer = ProteaseAnalyzer("protease_database.tsv")

# Analyze a single sequence
sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUSFLFIRKASMKAYHLGFQAACMKYLSMDTDQQQMRQRQKRYQNWR"

results = analyzer.analyze_sequence(sequence, "Example_Protein")
analyzer.print_results(results)

# Analyze a FASTA file
results = analyzer.analyze_fasta("proteins.fasta")
analyzer.print_results(results)
analyzer.export_results(results, "cleavage_results.csv")

# The analyzer automatically handles conversion from three-letter codes in the database
# to single-letter codes for sequence matching. For example:
# Database: "Phe" -> "F", "Tyr" -> "Y", "Ala" -> "A"
# Special cases like "TyI" (iodinated tyrosine) are handled as "Y"
"""