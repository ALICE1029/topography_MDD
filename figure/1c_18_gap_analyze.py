import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import cm, colors
import matplotlib as mpl
# 读取 .mat 文件
mat = scipy.io.loadmat(r'D:\study\sub2\liu\network_gap.mat')
diff = scipy.io.loadmat(r'D:\study\sub2\liu\average_diff.mat')

hc_network_gap = mat['hc_network_gap_new']
sub1_network_gap = mat['sub1_network_gap_new']
sub2_network_gap = mat['sub2_network_gap_new']
average_diff = diff['average_diff']  # shape: 18 x 2

# 动态颜色映射设置
# Step 1: 设置 colormap（对称色图）
cmap = mpl.cm.coolwarm  # 或者你自己定义的红-灰-蓝色图
vmax = np.max(np.abs(diff['average_diff']))
norm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)

# Step 2: 生成 sub1 / sub2 的颜色（根据 diff）
colors_sub1 = [cmap(norm(val)) for val in diff['average_diff'][:, 0]]
colors_sub2 = [cmap(norm(val)) for val in diff['average_diff'][:, 1]]

# 指定不显著的网络索引（灰色显示）
sub1_grey_idx = [1, 3, 5,6,8]
sub2_grey_idx = [5, 7]

# 创建图形和子图
fig, axs = plt.subplots(1, 1, figsize=(15, 15))

def plot_kde(ax, data, color, mean_line, xlim=(-40, 40)):
  #  sns.kdeplot(data, ax=ax, color=color, fill=True, legend=False, alpha=1)
    ax.axvline(mean_line, color='darkgrey')
    ax.set_yticks([])
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlim(xlim)

for i in range(18):
    # HC 分布图
  #  plot_kde(axs[0, i], hc_network_gap[:, i], 'lightgrey', np.mean(hc_network_gap[:, i]))



    # Sub2 分布图
    if i in sub2_grey_idx:
        color2 = 'lightgrey'
    else:
        color2 = cmap(norm(average_diff[i, 1]))  # sub2 vs HC

  #  plot_kde(axs[1, i], sub2_network_gap[:, i], colors_sub2[i], np.mean(hc_network_gap[:, i]))
    # Sub1 分布图
    if i in sub1_grey_idx:
        color1 = 'lightgrey'
    else:
        color1 = cmap(norm(average_diff[i, 0])) # sub1 vs HC

   # plot_kde(axs[2, i], sub1_network_gap[:, i], colors_sub1[i], np.mean(hc_network_gap[:, i]))
fig.subplots_adjust(bottom=0.15)  # 留出空间

sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.02])  # 可调
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
cbar.set_label("Mean Difference Relative to HC")
plt.tight_layout()
plt.savefig('gap_dynamic_colormap.eps', format='eps', dpi=300)
plt.show()