import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# for j in [13]:
#     print(j)
#     data = np.loadtxt('D:/study/sub2/liu/treatment/'+str(j)+'_mfd_treatmentcsu.txt', delimiter=',')
#     baseline = data[:, 0]
#     followup = data[:, 1]
#     #group = np.array([1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1])
#     group = np.array([	2,2	,1,	1	,2,	1,	2,	2	,2	,1,	1,	2,	2,	1,	1,	1	,2	,2	,2	,2	,2	,2	,1	,1	,1	,1	,1	,1	,2	,1	,2	,1	,2	,2	,1	,1	,1	,1	,1	,1	,1	,1,	1])
#     # 创建 DataFrame
#     df = pd.DataFrame({
#         'group': group,
#         'pre_treatment': baseline,
#         'post_treatment': followup
#     })
#
#     # 重组数据以适应小提琴图
#     df_melted = pd.melt(df, id_vars=['group'], value_vars=['pre_treatment', 'post_treatment'],
#                         var_name='Treatment', value_name='Value')
#
#     # 创建新的列以表示横轴位置
#     df_melted['Position'] = df_melted.apply(
#         lambda row: (0 if row['group'] == 1 else 2) if row['Treatment'] == 'pre_treatment' else (1 if row['group'] == 1 else 3),
#         axis=1
#     )
#
#     # 设置颜色
#     palette = {2: 'lightblue', 1: 'lightcoral'}
#
#     # 创建画布
#     fig, ax = plt.subplots(figsize=(8, 6))
#
#     # 绘制小提琴图
#     sns.violinplot(x='Position', y='Value', hue='group', data=df_melted,
#                    palette=palette)
#     #sns.boxplot(x='Position', y='Value', data=df_melted, color='gray', width=0.2, boxprops={'zorder': 2})
#     # 添加连线
#     for i in range(len(df)):
#         if group[i] == 1:
#             plt.plot([-0.2, 0.8], [baseline[i], followup[i]], color='gray', alpha=0.5)
#             plt.scatter([-0.2, 0.8], [baseline[i], followup[i]], color='red', s=20)  # 添加点
#         else:
#             plt.plot([2.2, 3.2], [baseline[i], followup[i]], color='gray', alpha=0.5)
#             plt.scatter([2.2, 3.2], [baseline[i], followup[i]], color='blue', s=20)  # 添加点
#     mean_pre_1 = df[df['group'] == 1]['pre_treatment'].mean()
#     mean_post_1 = df[df['group'] == 1]['post_treatment'].mean()
#     mean_pre_2 = df[df['group'] == 2]['pre_treatment'].mean()
#     mean_post_2 = df[df['group'] == 2]['post_treatment'].mean()
#     plt.plot([-0.2, 0.8], [mean_pre_1, mean_post_1], color='red', linestyle='--', linewidth=2)
#     plt.plot([2.2, 3.2], [mean_pre_2, mean_post_2], color='blue', linestyle='--', linewidth=2)
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#
#     plt.tick_params(labelsize =20,direction = 'in')
#     plt.savefig(str(j)+'_gap_changecsu.png', dpi=600)
#     plt.legend('', frameon=False)
#     # 调整 x 轴刻度和标签
#     plt.xticks([0,1,2,3],[])
#     plt.xlabel('')
#     plt.ylabel('')
#.show()


data = np.loadtxt('D:/study/sub2/liu/treatment/treatment_mfd.txt')
baseline = data[:, 0]
followup = data[:, 1]
group = np.array([1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1])
#group = np.array([	2,2	,1,	1	,2,	1,	2,	2	,2	,1,	1,	2,	2,	1,	1,	1	,2	,2	,2	,2	,2	,2	,1	,1	,1	,1	,1	,1	,2	,1	,2	,1	,2	,2	,1	,1	,1	,1	,1	,1	,1	,1,	1])
# 创建 DataFrame
df = pd.DataFrame({
    'group': group,
    'pre_treatment': baseline,
    'post_treatment': followup
})

# 重组数据以适应小提琴图
df_melted = pd.melt(df, id_vars=['group'], value_vars=['pre_treatment', 'post_treatment'],
                    var_name='Treatment', value_name='Value')

# 创建新的列以表示横轴位置
df_melted['Position'] = df_melted.apply(
    lambda row: (0 if row['group'] == 1 else 2) if row['Treatment'] == 'pre_treatment' else (1 if row['group'] == 1 else 3),
    axis=1
)

# 设置颜色
palette = {2: 'lightblue', 1: 'lightcoral'}

# 创建画布
fig, ax = plt.subplots(figsize=(8, 6))

# 绘制小提琴图
sns.violinplot(x='Position', y='Value', hue='group', data=df_melted,
               palette=palette)
#sns.boxplot(x='Position', y='Value', data=df_melted, color='gray', width=0.2, boxprops={'zorder': 2})
# 添加连线
for i in range(len(df)):
    if group[i] == 1:
        plt.plot([-0.2, 0.8], [baseline[i], followup[i]], color='gray', alpha=0.5)
        plt.scatter([-0.2, 0.8], [baseline[i], followup[i]], color='red', s=20)  # 添加点
    else:
        plt.plot([2.2, 3.2], [baseline[i], followup[i]], color='gray', alpha=0.5)
        plt.scatter([2.2, 3.2], [baseline[i], followup[i]], color='blue', s=20)  # 添加点
mean_pre_1 = df[df['group'] == 1]['pre_treatment'].mean()
mean_post_1 = df[df['group'] == 1]['post_treatment'].mean()
mean_pre_2 = df[df['group'] == 2]['pre_treatment'].mean()
mean_post_2 = df[df['group'] == 2]['post_treatment'].mean()
plt.plot([-0.2, 0.8], [mean_pre_1, mean_post_1], color='red', linestyle='--', linewidth=2)
plt.plot([2.2, 3.2], [mean_pre_2, mean_post_2], color='blue', linestyle='--', linewidth=2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend('', frameon=False)
plt.tick_params(labelsize =20,direction = 'in')
plt.savefig('gap_change_mfd.png', dpi=600)
plt.legend('', frameon=False)
# 调整 x 轴刻度和标签
plt.xticks([0,1,2,3],[])
plt.xlabel('')
plt.ylabel('')
