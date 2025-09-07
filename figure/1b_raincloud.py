
import h5py
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ptitprince as pt

# 读取MAT文件
file_data1 = h5py.File(r'D:\study\sub2\liu\adjusted_all_gap_remove65sex_new.mat')
file_data2 = h5py.File(r'D:\study\sub2\liu\adjusted_mdd_all_gap_remove65sex.mat')
data1 = np.array(file_data1['/ad_gap_new']).flatten()
data2 = np.array(file_data2['/adjusted_mdd_all_gap']).flatten()


# 创建数据框用于绘图
import pandas as pd

df = pd.DataFrame({
    'Group': ['adjusted_all_gap'] * len(data1) + ['adjusted_mdd_all_gap'] * len(data2),
    'Value': np.concatenate((data1, data2))
})

# 绘制雨云图
plt.figure(figsize=(8, 6))
ax = pt.RainCloud(x='Group', y='Value', data=df, palette=["#B0B0B0", "#FFC0CB"], bw=0.2, width_viol=0.4, ax=None, alpha=1)
ax.tick_params(labelsize = 15,direction = 'in')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(labelsize=25, direction='in')
plt.title('Raincloud Plot of Two Datasets')
plt.yticks(np.arange(-20,20, step=10))
plt.show()
