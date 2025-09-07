import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 读取 Excel 文件
file_path = r"D:\study\sub2\gene\sub1\sub1_posgene.xlsx"
# 读取 Excel 文件，将其存储在 DataFrame 中
df = pd.read_excel(file_path)

# 提取所需列的数据
count = df['Count']
percent = df['%']
log10_p = df['Log10(P)']

neg_log10_p = -log10_p

# 将 % 转化为 log2(x)，添加一个小的正数偏移量以避免 log2 计算错误
epsilon = 1e-10
log2_percent = np.log2(percent + epsilon)

# 对 -log10(p) 进行归一化处理，使其透明度在 0 到 1 之间
alpha = (neg_log10_p - np.min(neg_log10_p)) / (np.max(neg_log10_p) - np.min(neg_log10_p))

# 增大散点大小，这里将 count 乘以一个倍数，例如乘以 50
scatter_size = count * 30
plt.figure(figsize=(8, 10))
# 绘制火山图

# 绘制散点图
# c='darkred' 表示点的颜色为深红色
# s=scatter_size 表示点的大小由增大后的 count 列的值决定
# alpha=alpha 表示点的透明度由归一化后的 -log10(p) 列的值决定
# edgecolors='none' 表示不显示点的边缘颜色
sc = plt.scatter(log2_percent, neg_log10_p, c='darkblue', s=scatter_size, alpha=alpha, edgecolors='none')

# 添加标题和标签
plt.title('Volcano Plot', fontsize=16)
plt.xlabel('log2(Percentage)', fontsize=14)
plt.ylabel('-log10(p)', fontsize=14)

# 设置横轴范围
plt.xlim(0.5, 3)

# 不显示上边界线和右边界线
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 设置横轴和纵轴刻度字体为 Arial 14 号
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontname('Arial')
    label.set_fontsize(30)

# 显示图形
plt.show()