#!/usr/bin/env python3
"""
tRNA Fragment Analyzer

This tool identifies tRNA fragments by matching them against a reference database
of full tRNA sequences and provides comprehensive analysis capabilities.

Features:
- Sequence-based matching of fragments to full tRNAs
- Mass-based validation of matches where available
- Identification of fragment positions within full tRNAs
- Support for modified nucleotides in tRNA sequences
- Results organization by tRNA type, amino acid, and structural region
"""

import os
import re
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any, Union

# Define nucleotide mass constants
NUCLEOTIDE_MASSES = {
    # Standard nucleotides
    'A': 347.2212,
    'U': 324.1813,
    'G': 363.2206,
    'C': 323.1965,

    # Modified nucleotides commonly found in tRNA
    'I': 348.2233,  # Inosine
    'Ψ': 324.1813,  # Pseudouridine (same mass as U)
    'm1A': 361.2478,  # 1-methyladenosine
    'm2G': 377.2472,  # N2-methylguanosine
    'D': 325.1773,  # Dihydrouridine
    'R': 347.2212,  # Purine (A or G) - using A mass as approximation
    'Y': 324.1813,  # Pyrimidine (C, U) - using U mass as approximation
    'N': 339.6799,  # Any nucleotide (average mass)
    'J': 361.2478,  # Modified A
    'L': 347.2212,  # Modified A
    'O': 363.2206,  # Modified G
    'P': 363.2206,  # Modified G
    '7': 377.2472,  # Modified G
    '#': 363.2206,  # Modified G
    '"': 323.1965,  # Modified C
    'K': 324.1813,  # Modified U
    'B': 325.1773,  # Modified U
    'M': 323.1965,  # Modified C
    'H': 325.1773,  # Modified U
}

# Average mass for unknown modifications
AVERAGE_NUCLEOTIDE_MASS = 339.6799
WATER_MASS = 18.0106


def calculate_rna_mass(sequence: str, include_terminal: str = 'CCA') -> float:
    """
    Calculate theoretical mass of RNA sequence

    Args:
        sequence: RNA sequence
        include_terminal: Whether to include terminal CCA (CCA, CC, C, or none)

    Returns:
        Calculated mass in Daltons
    """
    # Filter out non-nucleotide characters and calculate mass
    total_mass = 0.0
    valid_nucleotides = 0

    for char in sequence:
        if char in NUCLEOTIDE_MASSES:
            total_mass += NUCLEOTIDE_MASSES[char]
            valid_nucleotides += 1
        elif re.match(r'[a-zA-Z]', char):
            # For unknown nucleotide codes, use average mass
            total_mass += AVERAGE_NUCLEOTIDE_MASS
            valid_nucleotides += 1

    # Add mass of terminal CCA if specified
    if include_terminal == 'CCA':
        total_mass += NUCLEOTIDE_MASSES['C'] * 2 + NUCLEOTIDE_MASSES['A']
        valid_nucleotides += 3
    elif include_terminal == 'CC':
        total_mass += NUCLEOTIDE_MASSES['C'] * 2
        valid_nucleotides += 2
    elif include_terminal == 'C':
        total_mass += NUCLEOTIDE_MASSES['C']
        valid_nucleotides += 1

    # Subtract water lost in phosphodiester bonds
    if valid_nucleotides > 1:
        total_mass -= WATER_MASS * (valid_nucleotides - 1)

    return total_mass


def determine_structural_region(start: int, length: int, full_length: int) -> str:
    """
    Determine the structural region of a tRNA fragment based on its position

    Args:
        start: Start position in the full tRNA
        length: Length of the fragment
        full_length: Length of the full tRNA

    Returns:
        Structural region name
    """
    # Approximate positions of tRNA structural elements as percentages
    # These are based on the cloverleaf structure of tRNA
    positions = {
        'acceptor_stem': {'start': 0, 'end': 0.1, 'end_part': {'start': 0.9, 'end': 1.0}},
        'd_arm': {'start': 0.1, 'end': 0.25},
        'anticodon_arm': {'start': 0.25, 'end': 0.45},
        'variable_loop': {'start': 0.45, 'end': 0.55},
        't_arm': {'start': 0.55, 'end': 0.8},
        'discriminator': {'start': 0.8, 'end': 0.9}
    }

    # Calculate relative positions
    start_percent = start / full_length
    end_percent = (start + length) / full_length

    # Determine region
    if (start_percent <= positions['acceptor_stem']['end'] or
            start_percent >= positions['acceptor_stem']['end_part']['start']):
        return "Acceptor Stem"
    elif (start_percent >= positions['d_arm']['start'] and
          end_percent <= positions['d_arm']['end']):
        return "D Arm"
    elif (start_percent >= positions['anticodon_arm']['start'] and
          end_percent <= positions['anticodon_arm']['end']):
        return "Anticodon Arm"
    elif (start_percent >= positions['variable_loop']['start'] and
          end_percent <= positions['variable_loop']['end']):
        return "Variable Loop"
    elif (start_percent >= positions['t_arm']['start'] and
          end_percent <= positions['t_arm']['end']):
        return "T Arm"
    elif (start_percent >= positions['discriminator']['start'] and
          end_percent <= positions['discriminator']['end']):
        return "Discriminator"
    else:
        return "Multiple Regions"


class TRNAFragmentAnalyzer:
    """Analyzer for tRNA fragments to identify their origins in reference tRNAs"""

    def __init__(self, mass_tolerance: float = 50):
        """
        Initialize the analyzer

        Args:
            mass_tolerance: Mass matching tolerance in PPM
        """
        self.mass_tolerance = mass_tolerance
        self.results = None

    def analyze(self, experimental_file: str, reference_file: str,
                reference_sheet: str = 'Sheet4') -> Dict:
        """
        Analyze tRNA fragments against reference database

        Args:
            experimental_file: Path to experimental data file
            reference_file: Path to reference database file
            reference_sheet: Sheet name in reference file

        Returns:
            Analysis results
        """
        print(f"Loading reference database from {reference_file}...")
        ref_df = pd.read_excel(reference_file, sheet_name=reference_sheet)
        print(f"Loaded {len(ref_df)} reference tRNA entries")

        print(f"Loading experimental data from {experimental_file}...")
        exp_df = pd.read_excel(experimental_file)
        print(f"Loaded {len(exp_df)} experimental tRNA fragments")

        # Match fragments to reference sequences
        results = self._match_fragments_to_reference(exp_df, ref_df)

        # Organize results by different categories
        organized_results = self._organize_results(results)

        self.results = {
            'summary': {
                'total_fragments': len(exp_df),
                'fragments_with_matches': results['fragments_with_matches'],
                'total_matches': results['total_matches'],
                'matches_by_type': results['matches_by_type'],
                'matches_by_amino_acid': results['matches_by_amino_acid']
            },
            'fragment_results': results['results'],
            'organized_results': organized_results
        }

        return self.results

    def _match_fragments_to_reference(self, exp_df: pd.DataFrame,
                                      ref_df: pd.DataFrame) -> Dict:
        """
        Match experimental fragments to reference tRNAs

        Args:
            exp_df: Experimental fragment data
            ref_df: Reference tRNA data

        Returns:
            Matching results
        """
        results = []
        match_count = 0
        fragments_with_matches = 0

        # Track matches by tRNA type and amino acid
        matches_by_type = {}
        matches_by_amino_acid = {}

        for _, exp_row in exp_df.iterrows():
            if not isinstance(exp_row.get('sequence'), str) or len(exp_row.get('sequence', '')) < 3:
                continue  # Skip very short sequences

            # Clean experimental sequence for matching
            clean_exp_seq = re.sub(r'[^A-Za-z]', '', exp_row['sequence']).upper()
            matches = []

            for _, ref_row in ref_df.iterrows():
                # Try both modified and unmodified sequences
                sequences = []
                if 'Sequence' in ref_row and isinstance(ref_row['Sequence'], str):
                    sequences.append({'type': 'modified', 'seq': ref_row['Sequence']})
                if 'Unmodified sequence' in ref_row and isinstance(ref_row['Unmodified sequence'], str):
                    sequences.append({'type': 'unmodified', 'seq': ref_row['Unmodified sequence']})

                for seq_entry in sequences:
                    full_seq = seq_entry['seq']
                    clean_ref_seq = re.sub(r'[^A-Za-z]', '', full_seq).upper()

                    # Check if experimental sequence is a substring of reference sequence
                    position = clean_ref_seq.find(clean_exp_seq)

                    if position >= 0:
                        # Calculate structural region of the match
                        structural_region = determine_structural_region(
                            position, len(clean_exp_seq), len(clean_ref_seq))

                        # Create the match entry
                        match = {
                            'tRNA_type': ref_row.get('tRNA type', "Unknown"),
                            'amino_acid': ref_row.get('AA', "Unknown"),
                            'name': ref_row.get('Name', "Unknown"),
                            'organism': ref_row.get('Organism', "Unknown"),
                            'full_sequence': full_seq,
                            'experimental_sequence': exp_row['sequence'],
                            'experimental_mass': exp_row.get('ending_mass'),
                            'position': position,
                            'sequence_type': seq_entry['type'],
                            'structural_region': structural_region,
                            'modified_nucleotides': "Yes" if ref_row.get('Modified') else "No"
                        }

                        # Mass validation if both masses are available
                        if exp_row.get('ending_mass') and ref_row.get('Mass_CCA'):
                            fragment_ratio = len(clean_exp_seq) / len(clean_ref_seq)
                            estimated_fragment_mass = ref_row['Mass_CCA'] * fragment_ratio
                            mass_difference = abs(exp_row['ending_mass'] - estimated_fragment_mass)
                            ppm_difference = (mass_difference / estimated_fragment_mass) * 1000000

                            match['estimated_mass'] = estimated_fragment_mass
                            match['mass_difference_ppm'] = ppm_difference
                            match['mass_validated'] = ppm_difference <= self.mass_tolerance

                        matches.append(match)

                        # Update statistics
                        trna_type = match['tRNA_type']
                        aa = match['amino_acid']

                        matches_by_type[trna_type] = matches_by_type.get(trna_type, 0) + 1
                        matches_by_amino_acid[aa] = matches_by_amino_acid.get(aa, 0) + 1

            if matches:
                fragments_with_matches += 1
                match_count += len(matches)

            results.append({
                'fragment': {
                    'sequence': exp_row['sequence'],
                    'mass': exp_row.get('ending_mass'),
                    'iteration': exp_row.get('n_iteration'),
                    'ladder': exp_row.get('ladder_number')
                },
                'matches': matches,
                'match_count': len(matches)
            })

        return {
            'results': results,
            'fragments_with_matches': fragments_with_matches,
            'total_matches': match_count,
            'matches_by_type': matches_by_type,
            'matches_by_amino_acid': matches_by_amino_acid
        }

    def _organize_results(self, results: Dict) -> Dict:
        """
        Organize results by tRNA type, amino acid, and structural region

        Args:
            results: Matching results

        Returns:
            Organized results
        """
        by_trna_type = {}
        by_amino_acid = {}
        by_structural_region = {}

        # Group fragments by different categories
        for result in results['results']:
            if not result['matches']:
                continue

            # Group by tRNA type
            for match in result['matches']:
                # By tRNA type
                trna_type = match['tRNA_type']
                if trna_type not in by_trna_type:
                    by_trna_type[trna_type] = []

                # Avoid duplicate fragments
                fragment_exists = any(
                    item['fragment']['sequence'] == result['fragment']['sequence'] and
                    item['fragment']['mass'] == result['fragment']['mass']
                    for item in by_trna_type[trna_type]
                )

                if not fragment_exists:
                    by_trna_type[trna_type].append({
                        'fragment': result['fragment'],
                        'match': match
                    })

                # By amino acid
                aa = match['amino_acid']
                if aa not in by_amino_acid:
                    by_amino_acid[aa] = []

                aa_fragment_exists = any(
                    item['fragment']['sequence'] == result['fragment']['sequence'] and
                    item['fragment']['mass'] == result['fragment']['mass']
                    for item in by_amino_acid[aa]
                )

                if not aa_fragment_exists:
                    by_amino_acid[aa].append({
                        'fragment': result['fragment'],
                        'match': match
                    })

                # By structural region
                region = match['structural_region']
                if region not in by_structural_region:
                    by_structural_region[region] = []

                region_fragment_exists = any(
                    item['fragment']['sequence'] == result['fragment']['sequence'] and
                    item['fragment']['mass'] == result['fragment']['mass']
                    for item in by_structural_region[region]
                )

                if not region_fragment_exists:
                    by_structural_region[region].append({
                        'fragment': result['fragment'],
                        'match': match
                    })

        return {
            'by_tRNA_type': by_trna_type,
            'by_amino_acid': by_amino_acid,
            'by_structural_region': by_structural_region
        }

    def export_results_to_excel(self, output_file: str) -> None:
        """
        Export results to an Excel file

        Args:
            output_file: Path to output Excel file
        """
        if not self.results:
            raise ValueError("No results available. Run analyze() first.")

        with pd.ExcelWriter(output_file) as writer:
            # Summary sheet
            summary_data = [
                ['tRNA Fragment Analysis Summary'],
                [''],
                ['Total fragments analyzed', self.results['summary']['total_fragments']],
                ['Fragments with matches', self.results['summary']['fragments_with_matches']],
                ['Total matches found', self.results['summary']['total_matches']],
                ['']
            ]

            summary_data.append(['Matches by tRNA Type'])
            summary_data.append(['tRNA Type', 'Count'])

            # Add tRNA type counts
            for trna_type, count in self.results['summary']['matches_by_type'].items():
                summary_data.append([trna_type, count])

            summary_data.append([''])
            summary_data.append(['Matches by Amino Acid'])
            summary_data.append(['Amino Acid', 'Count'])

            # Add amino acid counts
            for aa, count in self.results['summary']['matches_by_amino_acid'].items():
                summary_data.append([aa, count])

            summary_df = pd.DataFrame(summary_data)
            summary_df.to_excel(writer, sheet_name='Summary', header=False, index=False)

            # All fragments sheet
            fragments_data = []
            for result in self.results['fragment_results']:
                fragments_data.append([
                    result['fragment']['sequence'],
                    result['fragment'].get('mass'),
                    result['fragment'].get('iteration'),
                    result['fragment'].get('ladder'),
                    result['match_count']
                ])

            fragments_df = pd.DataFrame(
                fragments_data,
                columns=['Fragment Sequence', 'Experimental Mass',
                         'Iteration', 'Ladder', 'Match Count']
            )
            fragments_df.to_excel(writer, sheet_name='All Fragments', index=False)

            # Detailed matches sheet
            matches_data = []
            for result in self.results['fragment_results']:
                if not result['matches']:
                    continue

                for match in result['matches']:
                    matches_data.append([
                        result['fragment']['sequence'],
                        result['fragment'].get('mass'),
                        match['tRNA_type'],
                        match['amino_acid'],
                        match['name'],
                        match['organism'],
                        match['position'],
                        match['structural_region'],
                        match['sequence_type'],
                        match['modified_nucleotides']
                    ])

            matches_df = pd.DataFrame(
                matches_data,
                columns=[
                    'Fragment Sequence', 'Experimental Mass', 'tRNA Type', 'Amino Acid',
                    'tRNA Name', 'Organism', 'Position', 'Structural Region',
                    'Sequence Type', 'Modified Nucleotides'
                ]
            )
            matches_df.to_excel(writer, sheet_name='Detailed Matches', index=False)

            # By structural region sheet
            region_data = []
            for region, fragments in self.results['organized_results']['by_structural_region'].items():
                for item in fragments:
                    region_data.append([
                        region,
                        item['fragment']['sequence'],
                        item['fragment'].get('mass'),
                        item['match']['tRNA_type'],
                        item['match']['amino_acid']
                    ])

            region_df = pd.DataFrame(
                region_data,
                columns=['Structural Region', 'Fragment Sequence',
                         'Mass', 'tRNA Type', 'Amino Acid']
            )
            region_df.to_excel(writer, sheet_name='By Structural Region', index=False)

        print(f"Results exported to {output_file}")

    def get_summary_stats(self) -> Dict:
        """
        Get summary statistics of the analysis

        Returns:
            Dictionary with summary statistics
        """
        if not self.results:
            raise ValueError("No results available. Run analyze() first.")

        # Count fragments by structural region
        region_counts = {}
        for region, fragments in self.results['organized_results']['by_structural_region'].items():
            region_counts[region] = len(fragments)

        return {
            'total_fragments': self.results['summary']['total_fragments'],
            'fragments_with_matches': self.results['summary']['fragments_with_matches'],
            'percent_matched': (self.results['summary']['fragments_with_matches'] /
                                self.results['summary']['total_fragments'] * 100),
            'total_matches': self.results['summary']['total_matches'],
            'matches_by_type': self.results['summary']['matches_by_type'],
            'matches_by_amino_acid': self.results['summary']['matches_by_amino_acid'],
            'matches_by_structural_region': region_counts
        }


def main():
    """Main function to demonstrate tRNA fragment analysis"""
    # Example usage
    analyzer = TRNAFragmentAnalyzer(mass_tolerance=100)
    results = analyzer.analyze(
        experimental_file='tRNA_total_Deg_1_result.xlsx',
        reference_file='Copy of tRNA_DB_Combined_v3.xlsx',
        reference_sheet='Sheet2'
    )

    # Print summary statistics
    stats = analyzer.get_summary_stats()
    print("\n=== Analysis Summary ===")
    print(f"Total fragments analyzed: {stats['total_fragments']}")
    print(f"Fragments with matches: {stats['fragments_with_matches']} "
          f"({stats['percent_matched']:.1f}%)")
    print(f"Total matches found: {stats['total_matches']}")

    # Print tRNA type distribution
    print("\n=== Matches by tRNA Type ===")
    for trna_type, count in stats['matches_by_type'].items():
        print(f"{trna_type}: {count}")

    # Print amino acid distribution
    print("\n=== Matches by Amino Acid ===")
    for aa, count in stats['matches_by_amino_acid'].items():
        print(f"{aa}: {count}")

    # Print structural region distribution
    print("\n=== Matches by Structural Region ===")
    for region, count in stats['matches_by_structural_region'].items():
        print(f"{region}: {count}")

    # Export results to Excel
    analyzer.export_results_to_excel('tRNA_fragment_analysis_results.xlsx')

    # Show sample matches
    print("\n=== Sample Fragment Matches ===")
    sample_results = [r for r in results['fragment_results'] if r['matches']][:3]

    for result in sample_results:
        print(f"\nFragment: {result['fragment']['sequence']} "
              f"(Mass: {result['fragment'].get('mass', 'N/A')})")
        print(f"Found {result['match_count']} potential matches:")

        for i, match in enumerate(result['matches'][:3]):
            print(f"  {i + 1}. {match['tRNA_type']} ({match['amino_acid']}) - {match['organism']}")
            print(f"     Position: {match['position']}, Region: {match['structural_region']}")
            print(f"     tRNA: {match['name']}")

        if result['match_count'] > 3:
            print(f"     ... and {result['match_count'] - 3} more matches")


if __name__ == "__main__":
    main()