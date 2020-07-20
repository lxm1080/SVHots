# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 10:05:14 2019
注释enhancer
@author: Lv
"""
import indextable
#以上是主函数里的

def setEnhancer(start_list,hotspot_id_index, statistic_table, cancer,hotspotLenth):
    
    enhancerfilelist,type=indextable.indextable(cancer,'enhancer')
    
    all=start_list
    title=all[0]
    del all[0]
    
    lencol=len(all[0])
    for item in all: 
        item.append('')
        item.append(0)
        item.append('')
        item.append(0)

    enhancer_dict={}
    for item in all:
        item_name=item[hotspot_id_index]
        enhancer_dict[item_name]=[]
        
        
   
    enhancertitle=[]
    
    #注释本癌症相关的enhancer
    sum_enhancer_genome=0
    sum_enhancer_hotspot=0
    if len(enhancerfilelist)!=0 :     
        for i in range(len(enhancerfilelist)):
            #file=open('data/enhancer/stomach_differentially_expressed_enhancers.bed','r')
            file=open('data/enhancer/'+enhancerfilelist[i]+'_differentially_expressed_enhancers.bed','r')
            for line in file:
                sum_enhancer_genome+=1
                #line_split=re.split(r'[;,\s]\s*', line) 
                line_split=line.split('\t') 
                chr=line_split[0][3:]
                start=int(line_split[1])
                end=int(line_split[2])
                enhancer=line_split[3]
                #superenhancer=line_split[3]+line_split[4]      
            
                for item in all:
                    item_name=item[hotspot_id_index]
                    if (chr==item[0]) and ((start>item[1] and start<item[2]) or (end>item[1] and end<item[2])):
                        item[lencol]=item[lencol]+enhancer+', '
                        item[lencol+1]=item[lencol+1]+1
                        sum_enhancer_hotspot+=1   
                        enhancer_dict[item_name].append(start)
            enhancertitle=enhancertitle+[type+' related enhancer','num.'+type+' related enhancer']
            cellbase_list=[type+' related enhancer',sum_enhancer_hotspot,sum_enhancer_genome-sum_enhancer_hotspot]
            statistic_table.loc[len(statistic_table)]=cellbase_list
                    
                              
    #注释所有enhancer          
    sum_allenhancer_genome=0
    sum_allenhancer_hotspot=0
    file=open('data/all_enhancer_1.bed','r')
    #file=open('data/enhancer/kidney_differentially_expressed_enhancers.bed','r')
    #file=open('data/allSuperEnhancer.hg38.bed','r')
    for line in file:
        sum_allenhancer_genome+=1
        #line_split=re.split(r'[;,\s]\s*', line) 
        line_split=line.split('\t') 
        chr=line_split[0][3:]
        start=int(line_split[1])
        end=int(line_split[2])
        enhancer=line_split[3]
        #superenhancer=line_split[3]+line_split[4]      
    
        for item in all:
            
            if (chr==item[0]) and ((start>item[1] and start<item[2]) or (end>item[1] and end<item[2])):
                item[lencol+2]=item[lencol+2]+enhancer+', '
                item[lencol+3]=item[lencol+3]+1
                sum_allenhancer_hotspot+=1   
    file.close()       
                
    
    enhancertitle=enhancertitle+['all enhancer','num.all enhancer']
    title=title+enhancertitle
    all.insert(0,title)
    
    #cellbase_list=['enhancer',type+' related enhancer', 'other enhancer', sum_enhancer_hotspot, sum_enhancer_genome-sum_enhancer_hotspot,sum_allenhancer_hotspot,sum_allenhancer_genome-sum_allenhancer_hotspot]
    cellbase_list=['all enhancer',sum_allenhancer_hotspot,sum_allenhancer_genome-sum_allenhancer_hotspot]
    statistic_table.loc[len(statistic_table)]=cellbase_list
    return type+' related enhancer',enhancer_dict

'''
    ratio_enhancer_genome=(sum_enhancer_genome-sum_enhancer_hotspot)/(3000000000-hotspotLenth)*1000000
    ratio_enhancer_hotspot=sum_enhancer_hotspot/hotspotLenth*1000000
    print('胃癌相关的enhancer在热点区的密度是'+ str(ratio_enhancer_hotspot)+'/MB')
    print('胃癌相关的enhancer在背景区的密度是'+ str(ratio_enhancer_genome)+'/MB')
'''    
   
    