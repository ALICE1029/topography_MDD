
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from loess.loess_1d import loess_1d

sns.set(style='ticks')

# Define file paths and colors

files = [
    ("dendritedev_sub2_domi.xlsx", "dendritedev_sub2_nodomi.xlsx"),
    ("synapsedev_sub2_domi.xlsx", "synapsedev_sub2_nodomi.xlsx")
# ("myelination_sub1_domi.xlsx", "myelination_sub1_nodomi.xlsx"),
# ("axondev_sub1_domi.xlsx", "axondev_sub1_nodomi.xlsx")
]
colors = ["tab:red", "tab:blue"]
# Create subplots
plt.figure(figsize=(12, 6))

for  (domi_file, nodomi_file), color in zip( files, colors):
    # Read and plot domi data
    ax = plt.gca()
    df = pd.read_excel(f"D:/study/sub2/brain_age/brain_span/{domi_file}")
    age = df['age'].values.flatten()
    y = df['y'].values.flatten()
    xout, yout, weights = loess_1d(age, y, frac=0.5)
    sns.lineplot(x=xout, y=yout, linewidth=3, color=color, ci=None, ax=ax)

    # Read and plot nodomi data
    df = pd.read_excel(f"D:/study/sub2/brain_age/brain_span/{nodomi_file}")
    age = df['age'].values.flatten()
    y = df['y'].values.flatten()
    xout, yout, weights = loess_1d(age, y, frac=0.5)
    sns.lineplot(x=xout, y=yout, linewidth=3, color=color, linestyle='--', ci=None, ax=ax)

    # Customize each subplot
    ax.plot((5.8074, 5.8074), (0, 1), 'k--', linewidth=1)
    ax.plot((8.0553, 8.0553), (0, 1), 'k--', linewidth=1)
    ax.plot((9.3015, 9.3015), (0, 1), 'k--', linewidth=1)
    ax.plot((11.2621, 11.2621), (0, 1), 'k--', linewidth=1)
    ax.plot((12.3923, 12.3923), (0, 1), 'k--', linewidth=1)
    ax.set_xlim(5.7, 14)
    #ax.set_xticks([5.8074, 8.0553, 9.3015, 11.2621, 12.3923])
   # ax.set_xticklabels(["8W", "Birth", "1Y", "6Y", "14Y"])
    ax.set_ylim(-0.01, 1.01)
    ax.set_yticks([0, 0.5, 1])
    ax.spines['bottom'].set_linewidth(3)
    ax.spines['left'].set_linewidth(3)
    ax.spines['top'].set_linewidth(3)
    ax.spines['right'].set_linewidth(3)

    ax.tick_params(labelsize=30, direction='in')
    ax.set_yticklabels(["0", "0.5", "1"])
    plt.tight_layout()

    #plt.hold(True)
plt.show()
