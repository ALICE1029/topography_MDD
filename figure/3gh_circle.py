import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# 定义颜色数量
num_colors = 20

# 生成随机颜色列表
colors = list(mcolors.TABLEAU_COLORS.values())[:num_colors]

# 随机打乱颜色列表顺序
np.random.shuffle(colors)

# 根据需要重复颜色列表
colors_info = colors * (num_colors // len(colors) + 1)

# 取前 num_colors 个颜色作为最终颜色信息
final_colors = colors_info[:num_colors]
# 读取 Excel 文件
df = pd.read_excel('D:/study/sub2/gene/sub1/sub1_posgene.xlsx')

# 按照 p 值降序排序
df_sorted = df.sort_values(by='Log10(P)', ascending=False)


# 绘制泡泡图
fig = plt.figure(figsize=(8, 6))


plt.scatter(df_sorted['Log10(P)'], df_sorted['Description'], s=df_sorted['Count']*10, c=final_colors, alpha=0.5)
plt.xlabel('-Log10(P)')
plt.grid(True)
plt.tight_layout()  # 自动调整布局

plt.show()
