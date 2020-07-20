# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 10:05:14 2019

@author: Lv
"""

#以上是主函数里的
import re    
import indextable

def setSEA(input_list,hotspot_id_index, statistic_table, cancer):
    
    all=input_list
    #input_len=len(input_list[1])
    title=all[0]
    del all[0]
    
    lencol=len(all[0])
    for item in all: 
        item.append('')
        item.append(0)
        item.append('')
        item.append(0)
        
    
    
    sum_allSEA_genome=0
    sum_allSEA_hotspot=0
    
    
    #file=open('data/allSuperEnhancer.hg19.bed','r')
    file=open('data/all_superenhancer.bed','r')
    #file=open('data/allSuperEnhancer.hg38.bed','r')
    for line in file:
        sum_allSEA_genome+=1
        line_split=re.split(r'[;,\s]\s*', line) 
        chr=line_split[0][3:]
        start=int(line_split[1])
        end=int(line_split[2])
        superenhancer=line_split[3]+line_split[4]      
    
        for item in all:
            if (chr==item[0]) and ((start>item[1] and start<item[2]) or (end>item[1] and end<item[2])):
                item[lencol]=item[lencol]+superenhancer+', '
                item[lencol+1]=item[lencol+1]+1
                sum_allSEA_hotspot+=1
    enhancertitle=['all super enhancer','num.all super enhancer'] #标题
            
    cellbase_list=['all super enhancer',sum_allSEA_hotspot,sum_allSEA_genome-sum_allSEA_hotspot]
    statistic_table.loc[len(statistic_table)]=cellbase_list
    
    
    
    sum_SEA_genome=0
    sum_SEA_hotspot=0   
    superenhancer_dict={}        
    for item in all:
        item_name=item[hotspot_id_index]
        superenhancer_dict[item_name]=[]#superenhancer_dict初始化，热点名为键值，位置为value
    
   
    superenhancerfilelist,type=indextable.indextable(cancer,'superenhancer') #返回的第一个参数是列表，第二个是字符串
    
    if len(superenhancerfilelist)!=0 :     
        for i in range(len(superenhancerfilelist)):
            filepath=open('data/dbsuper/'+superenhancerfilelist[i]+'.bed','r')
        
     
            #stomochfile=open('data/dbsuper/Stomach Smooth Muscle.bed','r')
           
            for line in filepath:
                sum_SEA_genome+=1
                line_split=re.split(r'[;,\s]\s*', line) 
                chr=line_split[0][3:]
                start=int(line_split[1])
                end=int(line_split[2])
                superenhancer=line_split[3]+line_split[4]#+line_split[5]
            
                for item in all:#每一个热点区间
                    #item_name=item[input_len-1]
                    item_name=item[hotspot_id_index]
                   
                    if (chr==item[0]) and ((start>item[1] and start<item[2]) or (end>item[1] and end<item[2])):
                       
                        item[lencol+2]=item[lencol+2]+superenhancer+', '
                        item[lencol+3]=item[lencol+3]+1
                        sum_SEA_hotspot+=1
                        superenhancer_dict[item_name].append(start)
        enhancertitle=enhancertitle+[type+' related super enhancer','num.'+type+' related super enhancer']
        '''
        cellbase_list=['super_enhancer',type+' related super enhancer', 'other super_enhancer', sum_SEA_hotspot, sum_SEA_genome-sum_SEA_hotspot,sum_allSEA_hotspot,sum_allSEA_genome-sum_allSEA_hotspot]
        '''
        cellbase_list=[type+' related super enhancer',  sum_SEA_hotspot, sum_SEA_genome-sum_SEA_hotspot]
        statistic_table.loc[len(statistic_table)]=cellbase_list

       
    title=title+enhancertitle
    all.insert(0,title)     
    return type+' related super enhancer',superenhancer_dict#如果没有相关的superenhancer，返回的是空字典
    
    
   # return ratio_allSEA_genome,ratio_allSEA_hotspot,ratio_SEA_genome,ratio_SEA_hotspot,superenhancer_dict
'''             
        ratio_allSEA_genome=(sum_allSEA_genome-sum_allSEA_hotspot)/(3000000000-hotspotLenth)*1000000
        ratio_allSEA_hotspot=sum_allSEA_hotspot/hotspotLenth*1000000
        ratio_SEA_genome= (sum_SEA_genome-sum_SEA_hotspot)/(3000000000-hotspotLenth)*1000000
        ratio_SEA_hotspot=sum_SEA_hotspot/hotspotLenth*1000000    
        print('dbsuper里的所有超级增强子在热点区的密度是'+ str(ratio_allSEA_hotspot)+'/MB')
        print('dbsuper里的所有超级增强子在背景区的密度是'+ str(ratio_allSEA_genome)+'/MB')
        print('dbsuper里的指定超级增强子在热点区的密度是'+ str(ratio_SEA_hotspot)+'/MB')
        print('dbsuper里的指定超级增强子在背景区的密度是'+ str(ratio_SEA_genome)+'/MB')    
'''        


        
