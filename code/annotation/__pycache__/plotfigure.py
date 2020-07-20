# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 10:29:05 2019
先画density图，再画热点信息图，但是热点信息图还没更新，现在的代码放在untitled0
@author: Lv
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import matplotlib






#画ratio图
def plotratio(df,path):
    for i in range(len(df)):
        title=df['element name'][i]
        xlable= 'Pvalue='  +   str(df['pvalue(poisson test)'][i])  +  ','  +  'OR='   +   str(df['OR'][i])
        y=[ df['density in hotspot'][i], df['density in nonhotspot'][i] ]      
        plt.figure(figsize=(6,3))      
        rect1 = [0,0,1,1] # [左, 下, 宽, 高] 规定的矩形区域 （全部是0~1之间的数，表示比例）           
        ax1 = plt.axes(rect1)
        num_list1 = [y[0]]     # 纵坐标值1
        num_list2 = [y[1]]      # 纵坐标值2
        x = range(len(num_list1))      
        ax1.bar(x=x, height=num_list1, width=0.17, alpha=0.6, color='darkblue', label="hotspot")
        ax1.bar(x=[i + 0.25 for i in x], height=num_list2, width=0.17, alpha=0.7,color='lightblue', label="nonhotspot")
        #if isyname==1:
        #    axindex.set_ylabel("numbers/MB",fontsize=40)#纵坐标label
        
        #font_title = {'family' : 'Times New Roman','weight' : 'normal', 'size'   : 28,}
        ax1.set_title(title, {'family' : 'Times New Roman','weight' : 'normal', 'size': 33,}) #设置标题
        ax1.set_xlabel(xlable,{'family' : 'Times New Roman','weight' : 'normal', 'size': 23,})#横坐标label
        ax1.tick_params(axis='y', labelsize=20)#纵坐标的刻度字体大小
        ax1.set_xticks([])#横坐标刻度不显示
        ax1.spines['right'].set_color('none')
        ax1.spines['top'].set_color('none')
        ax1.spines['bottom'].set_color('none')#设置坐标轴颜色为无       
        ax1.set_ylabel("density per MB",{'family' : 'Times New Roman','weight' : 'bold', 'size': 26,})#纵坐标label
        leg = ax1.legend(fontsize=20,loc=[1,1])
        leg.get_frame().set_linewidth(0.0) #设置图例的边框消失
        plt.savefig(path+df['element name'][i]+'_ratio.jpg',dpi=200, bbox_inches='tight')
        #ax2.legend(fontsize=20,loc=[1,1]) #设置标注，并把标注放在右上角


#下面三个函数是画interval图
from matplotlib.ticker import FuncFormatter
#输入列表，输出列表中不为0的个数，不为0的百分比
def countpercent(inputlist):
    sumnum=len(inputlist)
    count=0
    for i in range(sumnum):
        if(inputlist[i]!=0):
            count+=1
    return [count, round(count/sumnum,2)]

def to_percent(temp, position):#转为百分比
  return '%1.0f'%(100*temp) + '%'

def plot_interval_percent(path,file,superenhancertitle,enhancertitle,drivergenetitle,eQTLtitle,gwastitle):
    #path='Results/lung/'
    #file='Results/lung/annotation_lung.csv'
    #superenhancertitle='lung related super enhancer'
    #enhancertitle='lung related enhancer'
    #drivergenetitle='LUAD driver gene'
    #eQTLtitle='lung related eQTL'
    #gwastitle='lung cancer snp'
    
    file=pd.read_csv(file)
    x=[] 
    y_percent=[] 
    y_count=[] 
    
    #gene
    genelist=file['num.gene']
    x.append('gene')
    y_count.append(countpercent(genelist)[0])
    y_percent.append(countpercent(genelist)[1])
    
    #cancer census gene
    CCGlist=file['num.cancer census gene']
    x.append('cancer census gene')
    y_count.append(countpercent(CCGlist)[0])
    y_percent.append(countpercent(CCGlist)[1])
    
    #all driver gene
    alldriverlist=file['all cancer driver gene']
    x.append('all cancer driver gene')
    y_count.append(countpercent(alldriverlist)[0])
    y_percent.append(countpercent(alldriverlist)[1])
    
    #driver gene
    driverlist=file['num.'+drivergenetitle]
    x.append(drivergenetitle)
    y_count.append(countpercent(driverlist)[0])
    y_percent.append(countpercent(driverlist)[1])
    
    #all SEA
    allSEAlist=file['num.all super enhancer']
    x.append('all super enhancer')
    y_count.append(countpercent(allSEAlist)[0])
    y_percent.append(countpercent(allSEAlist)[1])
    
    #SEA
    if  superenhancertitle in file.columns:
        SEAlist=file['num.'+superenhancertitle]
        x.append(superenhancertitle)
        y_count.append(countpercent(SEAlist)[0])
        y_percent.append(countpercent(SEAlist)[1])
        
    #all enhancer
    allenhancerlist=file['num.all enhancer']
    x.append('all enhancer')
    y_count.append(countpercent(allenhancerlist)[0])
    y_percent.append(countpercent(allenhancerlist)[1])
    
    #enhancer
    if  enhancertitle in file.columns:
        enhancerlist=file['num.'+enhancertitle]
        x.append(enhancertitle)
        y_count.append(countpercent(enhancerlist)[0])
        y_percent.append(countpercent(enhancerlist)[1])
        
        
    #all gwas snp
    allgwaslist=file['num.all gwas snp']
    x.append('all gwas snp')
    y_count.append(countpercent(allgwaslist)[0])
    y_percent.append(countpercent(allgwaslist)[1])
    
    #gwas snp
    if  gwastitle in file.columns:
        gwaslist=file['num.'+gwastitle]
        x.append(gwastitle)
        y_count.append(countpercent(gwaslist)[0])
        y_percent.append(countpercent(gwaslist)[1])
        
    #eQTL
    if  eQTLtitle in file.columns:
        eQTLtlist=file['num.'+eQTLtitle]
        x.append(eQTLtitle)
        y_count.append(countpercent(eQTLtlist)[0])
        y_percent.append(countpercent(eQTLtlist)[1])
     
    #title='distribution of element in hotspot'
    plt.figure(figsize=(4,6))
    rect1 = [0,0,1,1] # [左, 下, 宽, 高] 规定的矩形区域 （全部是0~1之间的数，表示比例）
    ax = plt.axes(rect1)
    barh=ax.barh(x, y_percent,alpha=0.5,color='black',height=0.75)
    #ax.set_title(title, {'family' : 'Times New Roman','weight' : 'normal', 'size': 33,})
       
    ax.tick_params(axis='y',labelsize=20)#纵坐标的刻度字体大小
    ax.tick_params(axis='x', labelsize=18)#纵坐标的刻度字体大小
    ax.set_xlabel('number of intervals in '+str(len(file))+' hotspot intervals' ,{'size': 20})#纵坐标label
    ax.spines['bottom'].set_color('none')#设置坐标轴颜色为无
    ax.spines['top'].set_color('none')#设置坐标轴颜色为无
    ax.spines['left'].set_color('none')#设置坐标轴颜色为无
    ax.spines['right'].set_color('none')#设置坐标轴颜色为无
        
    for b,ycount in zip(barh,y_count):
        width = b.get_width()
        ax.annotate(ycount,
        xy=(width, b.get_y() + b.get_height() / 2),
        xytext=(0, 3), # 3 points vertical offset
        textcoords="offset points",color='black',
        ha='left', va='center',fontsize=20)    
    ax.xaxis.set_major_formatter(FuncFormatter(to_percent))    
    plt.savefig(path+'intervalpercent.jpg',dpi=200, bbox_inches='tight')   




#画热图
def plot_distribution(path,file,superenhancertitle,enhancertitle,gwastitle, eQTLtitle,drivergenetitle):
    
    command='Rscript  code/annotation/drawhotspot.R'+' --path '+path+' --file '+file+' --superenhancertitle '+superenhancertitle.replace(' ','.')+' --enhancertitle '+enhancertitle.replace(' ','.')+' --gwastitle '+gwastitle.replace(' ','.')+' --eQTLtitle '+eQTLtitle.replace(' ','.')+' --drivergenetitle '+drivergenetitle.replace(' ','.')

    os.system(command)






#下面的函数是画OR对比图
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def plot_OR_comparison(typelist):
    data = pd.DataFrame(columns=['element name'])
    for i in range(len(typelist)):
        
        temp=pd.read_csv('Results/'+typelist[i]+'/cancer_specifity_statistic.csv')
        if i==0:
            data['element name']=temp['element name']
        data[typelist[i]]=temp['OR']
        
    #纵坐标
    element_name=[]
    for i in data['element name']:
        element_name.append(i)
    #横坐标
    typename=typelist
    
    OR=[]
    for i in range(len(data.index)):#data每一行
        temp=[]
        row=data.loc[i]
        for j in range(1,len(data.columns)):
            temp.append(round(row[j],2))
            #temp.append(row[j])
        OR.append(temp)
        
    OR=np.array(OR)
    
    fig, ax = plt.subplots(figsize=(7,6))
    
    im, cbar = heatmap(OR, element_name, typename, ax=ax,
                       cmap="YlGn", cbarlabel="OR")
    texts = annotate_heatmap(im)#, valfmt="{x:.1f} t")
    
    fig.tight_layout()
    plt.savefig('Results/ORcomparison.jpg',dpi=200, bbox_inches='tight')

'''
2019画
def drawsub(axindex,firstname,secondname,xname,isyname,islegend):
    num_list1 = [firstname]      # 纵坐标值1
    num_list2 = [secondname]      # 纵坐标值2
    x = range(len(num_list1))
    axindex.bar(x=x, height=num_list1, width=0.1, alpha=0.8, color='gold', label="hotspot")
    axindex.bar(x=[i + 0.15 for i in x], height=num_list2, width=0.1, color='indigo', label="nonhotspot")
    if isyname==1:
        axindex.set_ylabel("numbers/MB",fontsize=40)#纵坐标label
    axindex.set_xlabel(xname,fontsize=40)#横坐标label
    axindex.tick_params(axis='y', labelsize=40)#纵坐标的刻度字体大小
    axindex.set_xticks([])#横坐标刻度不显示
    if islegend==1:
        axindex.legend(fontsize=20) 
        
def plotratiobefore(name,ratio_allSEA_genome,ratio_allSEA_hotspot,ratio_SEA_genome,ratio_SEA_hotspot,ratio_allGene_hotspot,ratio_allGene_genome,ratio_alldrivergene_hotspot,ratio_alldrivergene_genome,ratio_drivergene_hotspot,ratio_drivergene_genome,ratio_cancercensusgene_hotspot,ratio_cancercensusgene_genome,ratio_GWAS_hotspot,ratio_GWAS_genome,ratio_eQTL_hotspot,ratio_eQTL_genome):


    plt.figure(figsize=(25,16))
    
    rect1 = [0,0.55,0.2,0.4] # [左, 下, 宽, 高] 规定的矩形区域 （全部是0~1之间的数，表示比例）
    rect2 = [0.25,0.55,0.2,0.4]
    rect3 = [0.5,0.55,0.2,0.4]
    rect4 = [0.75,0.55,0.2,0.4]
    rect5 = [0,0,0.2,0.4]
    rect6 = [0.25,0,0.2,0.4]
    rect7 = [0.5,0,0.2,0.4]
    rect8 = [0.75,0,0.2,0.4]
    
    rect1 = [0,0.6,0.15,0.35] # [左, 下, 宽, 高] 规定的矩形区域 （全部是0~1之间的数，表示比例）
    rect2 = [0.25,0.6,0.15,0.35]
    rect3 = [0.5,0.6,0.15,0.35]
    rect4 = [0.75,0.6,0.15,0.35]
    rect5 = [0,0.1,0.15,0.35]
    rect6 = [0.25,0.1,0.15,0.35]
    rect7 = [0.5,0.1,0.15,0.35]
    rect8 = [0.75,0.1,0.15,0.35]
    
    
    ax1 = plt.axes(rect1)
    ax2 = plt.axes(rect2)
    ax3 = plt.axes(rect3)
    ax4 = plt.axes(rect4)
    ax5 = plt.axes(rect5)
    ax6 = plt.axes(rect6)
    ax7 = plt.axes(rect7)
    ax8 = plt.axes(rect8)
    
    
    
    drawsub(ax1,ratio_allSEA_hotspot,ratio_allSEA_genome,'allSEA',1,0)
    drawsub(ax2,ratio_SEA_hotspot,ratio_SEA_genome,'SEA',0,0)
    drawsub(ax3,ratio_allGene_hotspot,ratio_allGene_genome,'Gene',0,0)
    drawsub(ax4,ratio_alldrivergene_hotspot,ratio_alldrivergene_genome,'alldrivergene',0,1)
    
    drawsub(ax5,ratio_drivergene_hotspot,ratio_drivergene_genome,'drivergene',1,0)
    drawsub(ax6,ratio_cancercensusgene_hotspot,ratio_cancercensusgene_genome,'cancercensusgene',0,0)
    drawsub(ax7,ratio_GWAS_hotspot,ratio_GWAS_genome,'GWAS',0,0)
    drawsub(ax8,ratio_eQTL_hotspot,ratio_eQTL_genome,'eQTL',0,0)
    
    plt.savefig('Results/'+name+'/ratio.jpg',dpi=200, bbox_inches='tight')
'''    

'''
2020.5.10
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plotratio(df,path):

    for i in range(len(df)):
    

        title1=df['special type'][i]
        title2='all'+' '+df['elementname'][i]
        xlable_1='Pvalue='+str(df['pvalue_special type'][i])+','+'OR='+str(df['OR of special type'][i])
        
        xlable_2='Pvalue='+str(df['pvalue_all type'][i])+','+'OR='+str(df['OR of all type'][i])
        
        
        
        y1=[ df['special type density in hotspot'][i], df['special type density in nonhotspot'][i] ]
        y2=[ df['all type density in hotspot'][i], df['all type density in nonhotspot'][i] ]
        
        
        
        
        
        plt.figure(figsize=(30,8))
        
        rect1 = [0,0.55,0.2,0.4] # [左, 下, 宽, 高] 规定的矩形区域 （全部是0~1之间的数，表示比例）
        rect2 = [0.25,0.55,0.2,0.4]
        
        ax1 = plt.axes(rect1)
        ax2 = plt.axes(rect2)
        
        
        num_list1 = [y1[0]]     # 纵坐标值1
        num_list2 = [y1[1]]      # 纵坐标值2
        x = range(len(num_list1))
        
        ax1.bar(x=x, height=num_list1, width=0.17, alpha=0.6, color='darkblue', label="hotspot")
        ax1.bar(x=[i + 0.25 for i in x], height=num_list2, width=0.17, alpha=0.7,color='lightblue', label="nonhotspot")
        #if isyname==1:
        #    axindex.set_ylabel("numbers/MB",fontsize=40)#纵坐标label
        
        #font_title = {'family' : 'Times New Roman','weight' : 'normal', 'size'   : 28,}
        ax1.set_title(title1, {'family' : 'Times New Roman','weight' : 'normal', 'size': 33,}) #设置标题
        ax1.set_xlabel(xlable_1,{'family' : 'Times New Roman','weight' : 'normal', 'size': 23,})#横坐标label
        ax1.tick_params(axis='y', labelsize=20)#纵坐标的刻度字体大小
        ax1.set_xticks([])#横坐标刻度不显示
        ax1.spines['right'].set_color('none')
        ax1.spines['top'].set_color('none')
        ax1.spines['bottom'].set_color('none')#设置坐标轴颜色为无
        
        ax1.set_ylabel("density per MB",{'family' : 'Times New Roman','weight' : 'bold', 'size': 26,})#纵坐标label
        #ax1.legend(loc=[1, 1])
        #if islegend==1:
        
        
        num_list1 = [y2[0]]     # 纵坐标值1
        num_list2 = [y2[1]]
        ax2.legend(fontsize=20) 
        
        ax2.bar(x=x, height=num_list1, width=0.17, alpha=0.6, color='darkblue', label="hotspot")
        ax2.bar(x=[i + 0.25 for i in x], height=num_list2, width=0.17,alpha=0.7, color='lightblue', label="nonhotspot")
        #if isyname==1:
        #    axindex.set_ylabel("numbers/MB",fontsize=40)#纵坐标label
        ax2.set_title(title2, {'family' : 'Times New Roman','weight' : 'normal', 'size': 33,}) #设置标题
        ax2.set_xlabel(xlable_2, {'family' : 'Times New Roman','weight' : 'normal', 'size': 23,})#横坐标label
        ax2.tick_params(axis='y', labelsize=20)#纵坐标的刻度字体大小
        ax2.set_xticks([])#横坐标刻度不显示
        #if islegend==1:
        
        ax2.spines['right'].set_color('none')
        ax2.spines['top'].set_color('none')
        ax2.spines['bottom'].set_color('none')#设置坐标轴颜色为无
        leg = ax2.legend(fontsize=20,loc=[1,1])
        leg.get_frame().set_linewidth(0.0) #设置图例的边框消失
        plt.savefig(path+df['elementname'][i]+'_ratio.jpg',dpi=200, bbox_inches='tight')
        #ax2.legend(fontsize=20,loc=[1,1]) #设置标注，并把标注放在右上角
'''

#传入画画
'''    
    #sample个数的图

    d = pd.read_csv('Results/rs3_2_annotation.csv', usecols=['no.samples','length.bp','hotspot.id', 'cancer census gene','number.bps'])

    plt.figure(figsize=(12,7))
    plt.xticks(rotation=90)
    plt.title('hotspots of RS3',fontsize=18)
    plt.ylabel('number of samples/breakpoints contributing to hotspot',fontsize=13)
    
    fig=plt.bar(d['hotspot.id'],d['number.bps'],color='indigo', label='breakpoints')    
    fig=plt.bar(d['hotspot.id'],d['no.samples'],color='gold', label='samples')    
    plt.legend(loc='upper right')
    plt.savefig("Results/no.sample.jpg",dpi=200, bbox_inches='tight')


#测试
    d = pd.read_csv('Results/rs3_2_annotation.csv', usecols=['no.samples','hotspot.id', 'cancer census gene','number.bps'])
    d=d.sort_values(['number.bps'],ascending=False)
    plt.figure(figsize=(12,7))
    plt.xticks(rotation=90)
    plt.title('hotspots of RS3',fontsize=18)
    plt.ylabel('number of samples/breakpoints contributing to hotspot',fontsize=13)
    
    fig=plt.bar(d['hotspot.id'],d['number.bps'],color='indigo', label='breakpoints')    
    fig=plt.bar(d['hotspot.id'],d['no.samples'],color='gold', label='samples')    
    plt.legend(loc='upper right')
    plt.savefig("Results/no.sample.jpg",dpi=200, bbox_inches='tight')
'''