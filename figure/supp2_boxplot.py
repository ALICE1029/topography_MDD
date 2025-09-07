import matplotlib.pyplot as plt
import numpy as np
import hdf5storage
import pandas as pd
import seaborn as sns

# Load the data
all_data = hdf5storage.loadmat(r'D:\study\sub2\draw\mae_site.mat')
data = all_data['mae_site_all']  # 替换为实际变量名plt.yticks(np.arange(np.min(data), np.max(data), step=0.5))

# 检查数据形状
print(data)  # 应该输出 (12, 10)
palette = sns.color_palette("coolwarm", n_colors=10)
# 绘制小提琴图（每一列是一个小提琴图）
plt.figure(figsize=(6, 6))
sns.boxplot(data=data, palette=palette, linewidth=1)  # 可以设置内核类型：'point'、'stick'等
sns.stripplot(data=data, color='grey', alpha=0.6, jitter=True, size=5)
sns.despine(top=True, right=True)

plt.yticks(np.arange(0,12, step=3))
plt.xticks(np.arange(1,10, step=1))
plt.xticks(fontsize=14, fontname='Arial')  # 设置横轴标签为 1-10，字体为 Arial，字号为 14

# 设置纵轴标签字体为 Arial 和字号为 14
plt.yticks(fontsize=14, fontname='Arial')
# 添加标签和标题
plt.xlabel('Networks')
plt.ylabel('Values')
plt.title('Violin Plot of 10 Networks')

# 显示图形
plt.show()