import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def create_trna_sequence_coverage_plot(sequences, coverage, title):
    """
    Create a tRNA coverage plot with nucleotide sequences

    Parameters:
    sequences: DataFrame with nucleotide sequences (A, U, G, C)
    coverage: DataFrame with boolean values for coverage (True where covered)
    title: str, title for the plot
    """
    fig, ax = plt.subplots(figsize=(20, 12))

    # Hide axes
    ax.set_xticks([])
    ax.set_yticks([])

    # Create a grid
    num_rows, num_cols = sequences.shape

    # Draw grid lines
    for i in range(num_rows + 1):
        ax.axhline(i, color='black', linewidth=0.5)
    for j in range(num_cols + 1):
        ax.axvline(j, color='black', linewidth=0.5)

    # Add sequence letters and colored backgrounds
    for i in range(num_rows):
        for j in range(num_cols):
            # Add nucleotide letter
            if pd.notna(sequences.iloc[i, j]):
                ax.text(j + 0.5, num_rows - i - 0.5, sequences.iloc[i, j],
                        ha='center', va='center', fontsize=8)

                # Add green background if covered
                if coverage.iloc[i, j]:
                    ax.add_patch(plt.Rectangle((j, num_rows - i - 1), 1, 1,
                                               facecolor='lightgreen', alpha=0.5))

    # Add column numbers on top
    for j in range(num_cols):
        ax.text(j + 0.5, num_rows, str(j + 1),
                ha='center', va='bottom', fontsize=8)

    # Add tRNA labels on left
    for i in range(num_rows):
        ax.text(-0.5, num_rows - i - 0.5, sequences.index[i],
                ha='right', va='center', fontsize=8)

    # Set title
    plt.title(title, pad=20)

    # Adjust layout
    plt.tight_layout()
    return plt


# Example usage with sample data
def generate_sample_data(num_trnas=30, num_positions=100):
    # Generate random sequences
    nucleotides = ['A', 'U', 'G', 'C']
    sequences = np.random.choice(nucleotides, size=(num_trnas, num_positions))
    coverage = np.random.choice([True, False], size=(num_trnas, num_positions), p=[0.7, 0.3])

    # Create DataFrames
    trna_names = [f'tRNA{i + 1}' for i in range(num_trnas)]
    positions = [str(i + 1) for i in range(num_positions)]

    seq_df = pd.DataFrame(sequences, index=trna_names, columns=positions)
    cov_df = pd.DataFrame(coverage, index=trna_names, columns=positions)

    return seq_df, cov_df


# Generate sample data
sequence_data, coverage_data = generate_sample_data()

# Create plot
create_trna_sequence_coverage_plot(
    sequence_data,
    coverage_data,
    'Coverage of cytoplasmic human tRNAs based on combined CCA, CC and C-tailed 3′ ladder fragments'
)

plt.show()