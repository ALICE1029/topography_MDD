import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os

# 可视化全局设置
plt.rcParams['font.family'] = "Arial"
plt.rcParams['font.size'] = 15  # 设置字体大小

# 输入和输出文件夹
input_folder = r'D:\study\sub2\scatterplot'
output_folder = r'D:\study\sub2\scatterplot\validation_loso'

# 确保输出文件夹存在
os.makedirs(output_folder, exist_ok=True)

# 遍历 1 到 10
for i in range(1, 11):
    for condition in ['pos_sub1center', 'neg_sub1center', 'pos_sub2center', 'neg_sub2center']:
        # 构造文件路径
        file_path = os.path.join(input_folder, f'{condition}{i}.txt')
        output_path = os.path.join(output_folder, f'{condition}{i}.png')

        # 读取数据
        try:
            file_data = pd.read_csv(file_path, sep='\t', header=None)
            file_data.columns = ["PLS weights: topography", "PLS weights: HDRS subscale"]

            # 绘制散点图
            sns.lmplot(
                x="PLS weights: topography",
                y="PLS weights: HDRS subscale",
                data=file_data,
                height=5, aspect=1.2,
                scatter_kws={"color": "#A9A9A9", "edgecolor": "darkgray"},
                line_kws={"color": "#BC3C28", "linewidth": 2}
            )

            # 获取当前轴对象
            ax = plt.gca()
            ax.set_facecolor("white")

            # 设置边框颜色
            for spine in ['top', 'bottom', 'left', 'right']:
                ax.spines[spine].set_color("black")

            # 调整刻度参数
            ax.tick_params(labelsize=15, direction='in')

            # 设置背景网格
            ax.grid(which='major', ls='--', c='grey', alpha=.3)

            # 去掉坐标轴标题
            ax.set_xlabel("")
            ax.set_ylabel("")

            # 隐藏右上边框
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            # 保存图像
            plt.savefig(output_path, format='png', dpi=300)
            plt.close()

            print(f"Saved: {output_path}")

        except Exception as e:
            print(f"Error processing {file_path}: {e}")
