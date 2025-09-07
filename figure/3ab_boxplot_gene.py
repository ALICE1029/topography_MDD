import matplotlib.pyplot as plt
import numpy as np
import hdf5storage
import pandas as pd
import seaborn as sns

# Load the data
all_data = hdf5storage.loadmat(r'D:\study\sub2\draw\permut_sub2.mat')
Myelin = [i for j in all_data['Myelin'] for i in j]
Axon = [i for j in all_data['Axon'] for i in j]
Dendrite = [i for j in all_data['Dendrite'] for i in j]
Synapse = [i for j in all_data['Synapse'] for i in j]

data = pd.DataFrame({
    "Myelin": Myelin,
    "Axon": Axon,
    "Dendrite": Dendrite,
    "Synapse": Synapse,
})

# Prepare the data for seaborn
data_melted = pd.melt(data, var_name='Structure', value_name='Expression')

# Plotting
fig, ax = plt.subplots(figsize=(6, 8))
sns.boxplot(x='Structure', y='Expression', data=data_melted, ax=ax,
            width=0.3, showfliers=False, palette=['lightgreen', 'thistle', 'mistyrose', 'lavender'],
            boxprops=dict(alpha=0.2), medianprops=dict(color='black', linewidth=5))

sns.stripplot(x='Structure', y='Expression', data=data_melted, ax=ax, jitter=True, palette=['lightgreen', 'thistle', 'mistyrose', 'lavender'], alpha=0.2, size=3)

# Specific y-values
#specific_y_values = [-1.57, -3.26, 1.54, 1.91]  # SUB1
specific_y_values = [-0.02, 0.56, 0.87, 1.40] #SUB2
colors = ['green', 'purple', 'red', 'blue']

# Scatter the specific y-values
for y_value, position, color in zip(specific_y_values, np.arange(len(data.columns)), colors):
    ax.scatter(position, y_value, color=color, s=150, zorder=3)

# Customize the plot
ax.tick_params(labelsize=20, direction='in')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
plt.ylabel('Expression difference')
plt.grid(linewidth=0.5)

plt.show()
