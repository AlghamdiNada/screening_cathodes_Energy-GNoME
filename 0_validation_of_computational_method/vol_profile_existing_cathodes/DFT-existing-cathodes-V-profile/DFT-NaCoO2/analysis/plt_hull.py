import pandas as pd
import matplotlib.pyplot as plt

# Load the data
df = pd.read_csv("hull.csv", skipinitialspace=True)

# Plot
plt.figure(figsize=(6, 4))
plt.scatter(df['x'], df['formation energy'], c='blue', alpha=0.7, edgecolors='k')

# Labels and title
plt.xlabel('x')
plt.ylabel('Formation Energy (eV)')
plt.title('Formation Energy vs. x')
plt.grid(True)
plt.tight_layout()

# Save or show
plt.savefig("formation_energy_plot.png", dpi=300)
plt.show()


