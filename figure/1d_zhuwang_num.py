import matplotlib.pyplot as plt
import numpy as np
import hdf5storage

# 导入数据
AtlasLoading_Mat = hdf5storage.loadmat(r'D:\study\sub2\liu\cres_sub2_2_remove65sexmask.mat')
num = (AtlasLoading_Mat['num']).transpose()

# 准备数据
#labels = ['DAN1', 'FPN2', 'SMN3', 'DMN4', 'DMN5', 'LIM6', 'SMN7', 'SUB8', 'DMN9', 'SMN10', 'VIS11', 'DAN12', 'FPN13', 'FPN14', 'DMN15', 'VAN16', 'DAN17', 'VIS18']
labels = ['DAN', 'FPN', 'SMN', 'DMN',  'LIM','SUB', 'VIS', 'VAN']
num_vars = len(labels)


# 根据第一列正值大小排序
sorted_indices = np.argsort(-num[:, 0])#CHANGE
sorted_labels = np.array(labels)[sorted_indices]
print(sorted_labels)
sorted_num1 = num[:, 0][sorted_indices]

sorted_num2 = num[:, 1][sorted_indices]
# 计算每个网络占总数的百分比
total_sum = np.sum(sorted_num1)  # 计算总数
sorted_num1 = (sorted_num1 / total_sum) * 100
total_sum = np.sum(sorted_num2)  # 计算总数
sorted_num2 = (sorted_num2 / total_sum) * 100
# sorted_labels = np.array(labels)
# sorted_num1 = num[:, 0]
# sorted_num2 = num[:, 1]
# 为了在蛛网图中闭合图形，需要在数据末尾再加上起点的数据
angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
angles += angles[:1]

y_data = list(sorted_num1) + [sorted_num1[0]]
y_data2 = list(sorted_num2) + [sorted_num2[0]]
sorted_labels = list(sorted_labels) + [sorted_labels[0]]

# 绘制蛛网图

fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))

ax.fill(angles, y_data, color='darkred', alpha=0.3, edgecolor='k')
#ax.fill(angles, y_data2, color='lightcoral', alpha=0.3, edgecolor='k')
#ax.fill(angles, y_data, color='blue', alpha=0.3, edgecolor='k')
#ax.fill(angles, y_data2, color='lightskyblue', alpha=0.3, edgecolor='k')
ax.set_yticks([max(y_data) * 0.5, max(y_data)])
print(max(y_data) * 0.5)
print(max(y_data))

print(max(y_data2) * 0.5)
print(max(y_data2))
ax.set_xticks(angles)
ax.set_rlabel_position(30)  # 调整角度位置，使标签不紧贴边界
ax.tick_params(labelsize=15,direction='in') # 0/15
#

ax.set_xticks(angles)
ax.set_xticklabels(sorted_labels)
ax.set_xticklabels([])
ax.set_yticklabels([])
plt.show()
