import matplotlib.pyplot as plt
import numpy as np
import hdf5storage
# 构建数据
from scipy import io

# #Import data
#AtlasLoading_Mat = hdf5storage.loadmat(r'D:\study\sub2\loadingweight.mat');

AtlasLoading_Mat = hdf5storage.loadmat(r'D:\study\sub2\liu\hamdsub2_weight_remove65\sum.mat');
num = (AtlasLoading_Mat['weight_sum'])
x_data = ['1', '2', '3', '4', '5', '6', '7', '8','9', '10', '11', '12', '13', '14', '15', '16', '17', '18']
y_data=list(range(18))
y_data2=list(range(18))
for i in range(0,18):
    y_data[i]=abs(num[i][0])
    y_data2[i] = abs(num[i][1]);
#AtlasLoading_Mat = hdf5storage.loadmat(r'D:\study\sub2\color.mat');
#gsr
#colors = ['red', 'forestgreen', 'm', 'm', 'royalblue','royalblue','hotpink','gold','royalblue','hotpink','slateblue','grey','gold','hotpink','forestgreen','red','hotpink','gold']
#ngsr
colors = ['LimeGreen', 'LemonChiffon', 'CornflowerBlue', 'pink', 'pink','Lavender','CornflowerBlue','lightgrey','pink','CornflowerBlue','Thistle','LimeGreen','LemonChiffon','LemonChiffon','pink','LightCoral','LimeGreen','Thistle']

colors=np.array(colors)
all=(np.array(y_data))+(np.array(y_data2));

arrIndex = np.array(all).argsort()
y_data = np.array(y_data)[arrIndex]
y_data2 = np.array(y_data2)[arrIndex]
x_data_new = np.array(x_data)[arrIndex]
colors_new = (colors)[arrIndex]
# 绘图
width = 0.8
fig, ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.bar(x_data_new, y_data, width, color=colors_new,hatch="\\\\",ec='k')
# 关键在bottom参数
plt.bar(x_data_new, y_data2, width, bottom=y_data, color=colors_new,ec='k')
ax.tick_params(labelsize=15, direction='in')
plt.xlabel("Network")
plt.ylabel("Sum of weights")
# 显示图例
#plt.legend()
plt.show()

