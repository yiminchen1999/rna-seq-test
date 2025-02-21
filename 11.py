import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Sample Data (You would replace this with your actual data loading code)
data = {
    'Position': np.random.randint(50, 70, 100),
    'Retention Time': np.linspace(10, 16, 100),
    'Mass': np.linspace(5750, 7800, 100),
    'Intensity': np.random.randint(2500, 100000, 100),
    'Nucleotide': np.random.choice(['A', 'C', 'G', 'U', 'm1A'], 100)
}

df = pd.DataFrame(data)

# Create a scatter plot
plt.figure(figsize=(10, 5))
scatter = plt.scatter(df['Mass'], df['Retention Time'], c=df['Intensity'], cmap='viridis', s=100, edgecolors='k', alpha=0.75)

# Add annotations
for i, txt in enumerate(df['Nucleotide']):
    plt.annotate(txt, (df['Mass'][i], df['Retention Time'][i]), textcoords="offset points", xytext=(0,10), ha='center')

# Adding a color bar
cbar = plt.colorbar(scatter)
cbar.set_label('Intensity')

# Labeling axes
plt.xlabel('Mass (Da)')
plt.ylabel('Retention Time (min)')
plt.title('Mass Spectrometry 3\'-Ladder Data')

plt.grid(True)
plt.show()
