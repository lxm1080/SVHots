# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 15:54:29 2020

@author: Lv
"""

import pandas as pd

import matplotlib.pyplot as plt

def plothot(name):
   
    df = pd.read_csv('Results/'+name+'/'+name+'_hotspot.csv')
    def my_density(length, number):
        return round(number/length*1000000)
    df['density'] = df.apply(lambda col: my_density(col['length.bp'], col['number.bps']), axis=1)
    df=df.sort_values(by="length.bp" , ascending=False) #按照热点长度排序
    
    #预处理
    xpre=df['hotspot.id'] 
    x=[]
    for item in xpre:
        templist=item.split('_')
        tempstr=templist[0]+':'+str(round(float(templist[1].split('M')[0])))+'MB'
        x.append(tempstr)
        
    y1=[]
    for item in df['length.bp']:
        y1.append(item/1000000)
    y2=df['no.samples'] 
    
    
    #开始画图
    title='hotspot information'
    plt.figure(figsize=(12,10))
    rect1 = [0,0.5,1,0.3] # [左, 下, 宽, 高] 规定的矩形区域 （全部是0~1之间的数，表示比例）
    rect2 = [0,0.06,1,0.3]
    ax1 = plt.axes(rect1)
    ax2 = plt.axes(rect2)
    
    
    rects1=ax1.bar(x,y1, width=0.8, alpha=0.9, color='#00007f', label="hotspot")
    ax2.bar(x,y2, width=0.8, alpha=0.9, color='#7f0000', label="hotspot")
    
    
    #font_title = {'family' : 'Times New Roman','weight' : 'normal', 'size' : 28,}
    #ax1.set_title(title, {'family' : 'Times New Roman','weight' : 'normal', 'size': 33,}) #设置标题
    
    #ax1.set_title(title, fontsize=20)
    ax1.set_ylabel("length of hotspot interval/MB",fontsize=20)#纵坐标label
    ax2.set_ylabel("number of samples",fontsize=20)#纵坐标label
    
    #ax1.set_xlabel(xlable,{'family' : 'Times New Roman','weight' : 'normal', 'size': 23,})#横坐标label
    #ax1.set_xlabel(xlable,fontsize=23)#横坐标label
    
    ax1.tick_params(axis='y', labelsize=20)#纵坐标的刻度字体大小
    ax2.tick_params(axis='y', labelsize=20)#纵坐标的刻度字体大小
    ax1.tick_params(axis='x', labelsize=12,rotation=90)#纵坐标的刻度字体大小
    ax2.set_xticks([])#横坐标刻度不显示
    ax1.spines['right'].set_color('none')
    ax1.spines['top'].set_color('none')
    ax1.spines['bottom'].set_color('none')#设置坐标轴颜色为无
    ax2.spines['right'].set_color('none')
    ax2.spines['top'].set_color('none')
    ax2.spines['bottom'].set_color('none')#设置坐标轴颜色为无
    #ax1.set_ylabel("density per MB",{'family' : 'Times New Roman','weight' : 'bold', 'size': 26,})#纵坐标label
    ax2.xaxis.set_ticks_position('top')   #将X坐标轴移到上面
    ax2.invert_yaxis()        
    
    def autolabel(rects,annot_text,ax):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect,text in zip(rects,annot_text):
            height = rect.get_height()
            ax.annotate(text,
                        xy=(rect.get_x()+rect.get_width() / 2 , height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', rotation=90,va='bottom')
    
    annot_text=df['density']
    autolabel(rects1,annot_text,ax1)
    
    plt.savefig('Results/'+name+'/hotspot_information.jpg',dpi=200, bbox_inches='tight')
    #plt.close()'''


