import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import codecs


# 可视化全局设置

plt.rcParams['font.family'] = "Arial"
plt.rcParams['font.size'] = 15# 设置字体大小为14

#可视化绘制

file_data = pd.read_csv(r'D:\study\sub2\scatterplot\pos_sub1split2.txt', sep='\t', header=None)
file_data.columns = ["PLS weights: topography","PLS weights:  HDRS subscale"]
#file_data.columns = ["Observed HDRS  change","Predicted HDRS  change"]
sns.lmplot(x="PLS weights: topography", y="PLS weights:  HDRS subscale", data=file_data,height=5, aspect=1.2,scatter_kws={"color":"#A9A9A9","edgecolor":"darkgray"}
    ,line_kws={"color":"#BC3C28","linewidth":2})
#sns.lmplot(x="Observed HDRS  change", y="Predicted HDRS  change", data=file_data,height=5, aspect=1.2,scatter_kws={"color":"#A9A9A9","edgecolor":"darkgray"}
     #,line_kws={"color":"#BC3C28","linewidth":2})
# 由于不支持ax参数，我们可以这样操作
ax = plt.gca()
# 定制化操作
ax.set_facecolor("white")
for spine in ['top','bottom','left','right']:
    ax.spines[spine].set_color("black")
ax.tick_params(labelsize = 15,direction = 'in')
ax.grid(which='major',ls='--',c='grey',alpha=.3)
# ax.set_xlim(left=-2,right=2)
#ax.set_ylim(bottom=-5,top=5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#titlefontdict = {"size":40,"color":"k",'family':'Arial'}
#ax.set_title('Example01 Of Linear Regression Scatter Plot ',titlefontdict,pad=15)
#text_font = {'family':'Arial','size':50,'weight':'bold','color':'black'}

plt.show()
