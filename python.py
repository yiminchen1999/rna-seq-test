import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_path = 'Phe_50FA_15pmol_T01_250114_LSS.xlsx'
data = pd.read_excel(file_path, sheet_name='Sheet1')

# Print column names to verify
print("Columns in the dataset:", data.columns)

# Preprocessing
# Use the correct column name for 'Monoisotopic Mass'
monoisotopic_mass_column = 'Monoisotopic Mass'  # Adjust this if the actual column name is different
if monoisotopic_mass_column not in data.columns:
    raise KeyError(f"Column '{monoisotopic_mass_column}' not found in the dataset.")

# Filter out invalid rows and reset the index
data_cleaned = data[data[monoisotopic_mass_column].notnull()].reset_index(drop=True)

# Ensure 'Sum Intensity' is numeric
data_cleaned['Sum Intensity'] = pd.to_numeric(data_cleaned['Sum Intensity'], errors='coerce')

# Drop rows with NaN in 'Sum Intensity' after conversion
data_cleaned = data_cleaned.dropna(subset=['Sum Intensity'])

# Step 1: Filter Top Peaks
sampling_num = 100  # Number of top peaks to sample
data_top = data_cleaned.nlargest(sampling_num, 'Sum Intensity')

# Step 2: Mass Sum Analysis
def mass_sum(data, full_mass, tolerance=0.1):
    """Identify mass pairs that sum to the given full mass."""
    results = []
    for i, mass1 in enumerate(data[monoisotopic_mass_column]):
        for j, mass2 in enumerate(data[monoisotopic_mass_column]):
            if i != j and abs(mass1 + mass2 - full_mass) <= tolerance:
                results.append((mass1, mass2))
    return pd.DataFrame(results, columns=['Mass1', 'Mass2'])

# Example: Perform Mass Sum analysis for a given full mass
full_mass = 24252.345  # Example full mass
data_mass_sum = mass_sum(data_top, full_mass)

# Step 3: Gap Filling
def gap_fill(data, mass_pairs, tolerance=0.5):
    """Fill gaps by finding masses near identified pairs."""
    gap_filled = []
    for _, row in mass_pairs.iterrows():
        for mass in data[monoisotopic_mass_column]:
            if abs(mass - row['Mass1']) <= tolerance or abs(mass - row['Mass2']) <= tolerance:
                gap_filled.append(mass)
    return pd.DataFrame(gap_filled, columns=['Filled Mass'])

data_gap_filled = gap_fill(data_top, data_mass_sum)

# Step 4: Visualization
def plot_mass_distribution(data, title):
    plt.figure(figsize=(10, 6))
    plt.hist(data[monoisotopic_mass_column], bins=30, alpha=0.7, label='Mass Distribution')
    plt.xlabel('Monoisotopic Mass')
    plt.ylabel('Frequency')
    plt.title(title)
    plt.legend()
    plt.show()

# Plot the original and gap-filled data
plot_mass_distribution(data_top, 'Top Mass Distribution')
plot_mass_distribution(data_gap_filled, 'Gap-Filled Mass Distribution')

# Save results
#data_mass_sum.to_csv('mass_sum_results.csv', index=False)
#data_gap_filled.to_csv('gap_filled_results.csv', index=False)

print("Analysis complete. Results saved as CSV files.")
