#!/usr/bin/env python3

"""
tRNA Fragment Analyzer - Usage Example

This script demonstrates how to use the tRNA Fragment Analyzer
to identify tRNA fragments from mass spectrometry data and perform
comprehensive analysis with visualizations.

Features:
1. tRNA fragment identification by sequence matching
2. Mass-based validation of matches
3. Positional and structural region analysis
4. Sequence composition analysis
5. Visualization of results with multiple plots
6. Export to Excel for further analysis
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from trna_fragment_analyzer import TRNAFragmentAnalyzer

# Set plotting style
plt.style.use('ggplot')
sns.set(style="whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12


def run_analysis():
    """Run tRNA fragment analysis and visualize results"""
    print("Starting tRNA fragment analysis...")

    # Initialize the analyzer
    analyzer = TRNAFragmentAnalyzer(mass_tolerance=50)

    # Run the analysis
    results = analyzer.analyze(
        experimental_file='tRNA_total_Deg_1_result.xlsx',
        reference_file='Copy of tRNA_DB_Combined_v3.xlsx',
        reference_sheet='Sheet2'
    )

    # Get summary statistics
    stats = analyzer.get_summary_stats()

    # Print summary statistics
    print("\n=== Analysis Summary ===")
    print(f"Total fragments analyzed: {stats['total_fragments']}")
    print(f"Fragments with matches: {stats['fragments_with_matches']} "
          f"({stats['percent_matched']:.1f}%)")
    print(f"Total matches found: {stats['total_matches']}")

    # Export results to Excel
    output_file = '/Users/chenyimin/PycharmProjects/rna-seq-test/tRNA_fragment_analysis_results.xlsx'
    analyzer.export_results_to_excel(output_file)
    print(f"\nResults exported to {output_file}")

    # Create visualizations
    create_visualizations(stats, results)

    # Additional analyses
    mass_stats = analyze_mass_distribution(results)
    position_data = analyze_fragment_positions(results)
    sequence_analysis = analyze_sequence_composition(results)

    return results


def create_visualizations(stats, results):
    """Create visualizations of analysis results"""
    print("\nGenerating visualizations...")
    os.makedirs('plots', exist_ok=True)

    # Figure 1: Fragment match distribution
    plt.figure(figsize=(10, 6))
    plt.bar(['Matched', 'Unmatched'],
            [stats['fragments_with_matches'],
             stats['total_fragments'] - stats['fragments_with_matches']])
    plt.title('Fragment Match Distribution')
    plt.ylabel('Number of Fragments')
    plt.tight_layout()
    #plt.savefig('plots/fragment_match_distribution.png')

    # Figure 2: Matches by tRNA type
    if stats['matches_by_type']:
        plt.figure(figsize=(12, 8))
        types_df = pd.DataFrame({
            'tRNA Type': list(stats['matches_by_type'].keys()),
            'Count': list(stats['matches_by_type'].values())
        })
        types_df = types_df.sort_values('Count', ascending=False)

        sns.barplot(x='Count', y='tRNA Type', data=types_df)
        plt.title('Matches by tRNA Type')
        plt.tight_layout()
        #plt.savefig('plots/matches_by_trna_type.png')

    # Figure 3: Matches by amino acid
    if stats['matches_by_amino_acid']:
        plt.figure(figsize=(12, 8))
        aa_df = pd.DataFrame({
            'Amino Acid': list(stats['matches_by_amino_acid'].keys()),
            'Count': list(stats['matches_by_amino_acid'].values())
        })
        aa_df = aa_df.sort_values('Count', ascending=False)

        sns.barplot(x='Count', y='Amino Acid', data=aa_df)
        plt.title('Matches by Amino Acid')
        plt.tight_layout()
        #plt.savefig('plots/matches_by_amino_acid.png')

    # Figure 4: Matches by structural region
    if stats['matches_by_structural_region']:
        plt.figure(figsize=(10, 6))
        region_df = pd.DataFrame({
            'Structural Region': list(stats['matches_by_structural_region'].keys()),
            'Count': list(stats['matches_by_structural_region'].values())
        })
        region_df = region_df.sort_values('Count', ascending=False)

        plt.pie(region_df['Count'], labels=region_df['Structural Region'],
                autopct='%1.1f%%', startangle=90)
        plt.axis('equal')
        plt.title('Fragment Distribution by Structural Region')
        plt.tight_layout()
        #plt.savefig('plots/matches_by_structural_region.png')

    # Figure 5: Fragment length distribution
    fragment_lengths = [len(result['fragment']['sequence'])
                        for result in results['fragment_results']]

    plt.figure(figsize=(10, 6))
    sns.histplot(fragment_lengths, bins=20, kde=True)
    plt.title('Fragment Length Distribution')
    plt.xlabel('Fragment Length (nucleotides)')
    plt.ylabel('Count')
    plt.tight_layout()
    #plt.savefig('plots/fragment_length_distribution.png')

    # Figure 6: Match counts per fragment
    match_counts = [result['match_count'] for result in results['fragment_results']]

    plt.figure(figsize=(10, 6))
    sns.histplot(match_counts, bins=20, kde=True)
    plt.title('Match Counts per Fragment')
    plt.xlabel('Number of Matches')
    plt.ylabel('Count')
    plt.tight_layout()
    #plt.savefig('plots/match_counts_distribution.png')

    print(f"Visualizations saved to the 'plots' directory")


def analyze_mass_distribution(results):
    """Analyze mass distribution of fragments"""
    print("\nAnalyzing mass distribution...")

    # Extract masses
    masses = []
    for result in results['fragment_results']:
        if result['fragment'].get('mass') is not None:
            masses.append(result['fragment']['mass'])

    if not masses:
        print("No mass data available")
        return

    # Calculate statistics
    mass_stats = {
        'count': len(masses),
        'min': min(masses),
        'max': max(masses),
        'mean': sum(masses) / len(masses),
        'bins': []
    }

    # Create mass bins
    bin_size = 1000  # 1 kDa bins
    bins = {}
    for mass in masses:
        bin_key = int(mass // bin_size) * bin_size
        if bin_key not in bins:
            bins[bin_key] = 0
        bins[bin_key] += 1

    mass_stats['bins'] = bins

    # Print statistics
    print(f"Mass range: {mass_stats['min']:.2f} - {mass_stats['max']:.2f} Da")
    print(f"Mean mass: {mass_stats['mean']:.2f} Da")

    # Plot mass distribution
    plt.figure(figsize=(12, 6))

    bin_keys = sorted(bins.keys())
    bin_values = [bins[key] for key in bin_keys]

    plt.bar([str(k) for k in bin_keys], bin_values)
    plt.title('Mass Distribution of tRNA Fragments')
    plt.xlabel('Mass (Da)')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.tight_layout()
    #plt.savefig('plots/mass_distribution.png')

    return mass_stats


def analyze_fragment_positions(results):
    """Analyze positions of fragments within full tRNA sequences"""
    print("\nAnalyzing fragment positions...")

    # Collect position data
    positions = []
    for result in results['fragment_results']:
        for match in result['matches']:
            if 'position' in match:
                positions.append({
                    'fragment': result['fragment']['sequence'],
                    'position': match['position'],
                    'tRNA_type': match['tRNA_type'],
                    'amino_acid': match['amino_acid'],
                    'region': match['structural_region']
                })

    if not positions:
        print("No position data available")
        return

    # Convert to DataFrame for analysis
    pos_df = pd.DataFrame(positions)

    # Create position heatmap by tRNA type and amino acid
    pivot_data = pos_df.pivot_table(
        index='tRNA_type',
        columns='region',
        values='position',
        aggfunc='count',
        fill_value=0
    )

    plt.figure(figsize=(14, 10))
    sns.heatmap(pivot_data, annot=True, cmap='YlGnBu', fmt='d')
    plt.title('Fragment Positions by tRNA Type and Structural Region')
    plt.tight_layout()
    plt.show()
    #plt.savefig('plots/position_heatmap.png')

    # Position histogram
    plt.figure(figsize=(12, 6))
    sns.histplot(data=pos_df, x='position', hue='region', multiple='stack', bins=20)
    plt.title('Distribution of Fragment Positions within tRNAs')
    plt.xlabel('Position in Full tRNA')
    plt.ylabel('Count')
    plt.tight_layout()
    #plt.savefig('plots/position_distribution.png')

    return pos_df


def analyze_sequence_composition(results):
    """Analyze nucleotide composition of tRNA fragments"""
    print("\nAnalyzing sequence composition...")

    # Collect sequence data
    sequences = []
    for result in results['fragment_results']:
        if result['matches']:
            frag_seq = result['fragment']['sequence']
            sequences.append({
                'sequence': frag_seq,
                'length': len(frag_seq),
                'has_match': True,
                'match_count': result['match_count']
            })
        else:
            frag_seq = result['fragment']['sequence']
            sequences.append({
                'sequence': frag_seq,
                'length': len(frag_seq),
                'has_match': False,
                'match_count': 0
            })

    if not sequences:
        print("No sequence data available")
        return

    # Convert to DataFrame
    seq_df = pd.DataFrame(sequences)

    # Calculate nucleotide frequencies
    nucleotides = ['A', 'U', 'G', 'C', 'Other']
    for nucleotide in nucleotides:
        if nucleotide != 'Other':
            seq_df[f'{nucleotide}_count'] = seq_df['sequence'].apply(
                lambda seq: seq.upper().count(nucleotide))
            seq_df[f'{nucleotide}_freq'] = seq_df[f'{nucleotide}_count'] / seq_df['length']
        else:
            # Count non-standard nucleotides
            seq_df['Other_count'] = seq_df.apply(
                lambda row: row['length'] - (row['A_count'] + row['U_count'] +
                                             row['G_count'] + row['C_count']),
                axis=1
            )
            seq_df['Other_freq'] = seq_df['Other_count'] / seq_df['length']

    # Calculate GC content
    seq_df['GC_content'] = (seq_df['G_count'] + seq_df['C_count']) / seq_df['length']

    # Nucleotide composition plot
    avg_freqs = [seq_df[f'{n}_freq'].mean() for n in nucleotides]

    plt.figure(figsize=(10, 6))
    plt.bar(nucleotides, avg_freqs)
    plt.title('Average Nucleotide Composition of tRNA Fragments')
    plt.xlabel('Nucleotide')
    plt.ylabel('Average Frequency')
    plt.ylim(0, 0.5)  # Adjust as needed
    for i, v in enumerate(avg_freqs):
        plt.text(i, v + 0.01, f'{v:.3f}', ha='center')
    plt.tight_layout()
    #plt.savefig('plots/nucleotide_composition.png')

    # GC content distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(data=seq_df, x='GC_content', hue='has_match')
    plt.title('GC Content Distribution of tRNA Fragments')
    plt.xlabel('GC Content')
    plt.ylabel('Count')
    plt.tight_layout()
    #plt.savefig('plots/gc_content_distribution.png')

    # Length vs. Match count scatter plot
    plt.figure(figsize=(10, 6))
    matched_df = seq_df[seq_df['has_match']].copy()
    plt.scatter(matched_df['length'], matched_df['match_count'], alpha=0.6)
    plt.title('Fragment Length vs. Number of Matches')
    plt.xlabel('Fragment Length (nucleotides)')
    plt.ylabel('Number of Matches')
    plt.tight_layout()
    #plt.savefig('plots/length_vs_matches.png')

    return seq_df

analyzer = TRNAFragmentAnalyzer(mass_tolerance=100)
results = analyzer.analyze(
    experimental_file='/Users/chenyimin/PycharmProjects/rna-seq-test/tRNA_total_Deg_1_result.xlsx',
    reference_file='/Users/chenyimin/PycharmProjects/rna-seq-test/Copy of tRNA_DB_Combined_v3.xlsx',
    reference_sheet='Sheet2'
)


stats = analyzer.get_summary_stats()
analyzer.export_results_to_excel('/Users/chenyimin/PycharmProjects/rna-seq-test/analysis_results.xlsx')
stats
results
analyze_mass_distribution(results)
analyze_sequence_composition(results)
create_visualizations(stats, results)
run_analysis()