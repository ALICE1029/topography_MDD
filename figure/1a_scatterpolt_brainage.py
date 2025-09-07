import matplotlib.pyplot as plt
import os
import hdf5storage
import numpy as np
import seaborn as sns
from scipy.stats import zscore
font = {'family': 'Arial', 'size': 14}

# 创建正方形画布
fig, ax = plt.subplots(figsize=(8, 8), dpi=600)
color = ['bisque', 'burlywood', 'orange', 'khaki', 'plum', 'thistle', 'violet', 'lavender', 'lightsteelblue', 'slategrey']
X = []
Y = []
group = np.repeat(np.arange(1, 11), [248, 108, 100, 150, 80, 200, 180, 300, 120, 314])

# 加载数据并创建散点图
for i in range(1, 11):
    FileName = f'mdd_{i}xy.mat'
    path = 'D:\\study\\sub2\\liu\\plot'
    File = os.path.join(path, FileName)
    AtlasLoading_Mat = hdf5storage.loadmat(File)

    x = AtlasLoading_Mat['x_tmp'].flatten()
    y = AtlasLoading_Mat['y_tmp'].flatten()
    X.extend(x)
    Y.extend(y)
    labelt = f'site{i}'

    ax.scatter(x, y, c=color[i - 1], s=80, alpha=0.4, label=labelt)

# 添加95%置信区间的拟合线

# X = zscore(X)
# Y = zscore(Y)
sns.regplot(x=X, y=Y, scatter=False, ax=ax, ci=95, color="black", line_kws={'linewidth': 2})


# 定制化图像外观
ax.set_aspect('equal', adjustable='box')  # 确保坐标轴比例相同，保持正方形
ax.set_facecolor("white")
for spine in ['top', 'bottom', 'left', 'right']:
    ax.spines[spine].set_color("black")
ax.tick_params(labelsize=30, direction='in')
ax.grid(which='major', ls='--', c='grey', alpha=.3)
ax.set_xlabel("", fontsize=30)
ax.set_ylabel("", fontsize=30)
# ax.set_xlim(left=8, right=65)  # Set x-axis limit
# ax.set_ylim(bottom=10, top=65)  # Set y-axis limit
# 设置横坐标范围，确保不超过数据的最小值和最大值
ax.set_xlim(5, 75)
ax.set_ylim(0, 75)  # 也可以对Y轴进行同样的操作

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 保存为PNG格式
plt.savefig("scatter_with_fit_line.png", format='png', bbox_inches='tight')
plt.show()
