# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 09:09:31 2019

@author: Lv
"""
import re    
import indextable

def setGene(input_list,hotspot_id_index,statistic_table,cancer,hotspotLenth):
    driverfile=indextable.indextable(cancer,'drivergene')
    
    all=input_list
    title=all[0]
    del all[0]
    
    lencol=len(all[0])
    for item in all: 
        item.append('')
        item.append(0)
        item.append('')
        item.append(0)
        item.append('')
        item.append(0)
        item.append('')
        item.append(0)
 
    GenePos = open('data/hg19_name.txt','r')

      
    #全部驱动基因字典    
    allDrivergene = open('data/allDrivergene.txt','r')
    allDrivergene_list=allDrivergene.read()
    allDrivergene_list = re.split(r'[;,\s]\s*', allDrivergene_list)    
    allDrivergene_dict={}
    for item in allDrivergene_list:
        if item not in allDrivergene_dict.keys():
            allDrivergene_dict[item]=1
        else:
            allDrivergene_dict[item]+=1
            
            
    #cancer cesus 基因列表
    CCG=open('data/CCG.txt','r')
    CCG_list=CCG.read()
    CCG_list = re.split(r'[;,\s]\s*', CCG_list)    
    

    
    #前列腺癌的字典
    Stomach_adenocarcinoma_DG=open('data/driverDB/'+driverfile+'(TCGA,US).txt','r')
    Stomach_adenocarcinoma_DG_list=Stomach_adenocarcinoma_DG.read()
    Stomach_adenocarcinoma_DG_list = re.split(r'[;,\s]\s*', Stomach_adenocarcinoma_DG_list)    
    Stomach_adenocarcinoma_DG_dict={}
    for item in Stomach_adenocarcinoma_DG_list:
        if item not in Stomach_adenocarcinoma_DG_dict.keys():
            Stomach_adenocarcinoma_DG_dict[item]=1
        else:
            Stomach_adenocarcinoma_DG_dict[item]+=1


    #开始注释
    sum_allGene_genome=0
    sum_alldrivergene_genome=len(allDrivergene_dict)
    sum_drivergene_genome=len(Stomach_adenocarcinoma_DG_dict)
    sum_cancercensusgene_genome=len(CCG_list)
    
    sum_allGene_hotspot=0
    sum_alldrivergene_hotspot=0
    sum_drivergene_hotspot=0
    sum_cancercensusgene_hotspot=0
    

    gene_dict={}
    keygene_dict={}
    for item in all:
           item_name=item[hotspot_id_index]
           gene_dict[item_name]=[]
           keygene_dict[item_name]=[]
           
    for line in GenePos:
        
        sum_allGene_genome+=1
        line_split=re.split(r'[;,\s]\s*', line) 
        chr=line_split[1][3:]
        start=int(line_split[2])
        end=int(line_split[3])
        gene=line_split[0]
    
        for item in all:
            
            item_name=item[hotspot_id_index]
            
            if (chr==item[0]) and ((start>item[1] and start<item[2]) or (end>item[1] and end<item[2])):
                item[lencol]=item[lencol]+gene+','
                item[lencol+1]=item[lencol+1]+1
                sum_allGene_hotspot+=1
                
                gene_dict[item_name].append(start)
     
                
                if gene in allDrivergene_dict.keys(): #and allDrivergene_dict[gene]>3:
                #if gene in allDrivergene:
                    item[lencol+2]=item[lencol+2]+gene+', '
                    item[lencol+3]=item[lencol+3]+1
                    sum_alldrivergene_hotspot+=1
                                     
                    
                if gene in Stomach_adenocarcinoma_DG_dict.keys():# and Prostate_adenocarcinoma_DG_dict[gene]>1:
                    item[lencol+4]=item[lencol+4]+gene+', '
                    item[lencol+5]=item[lencol+5]+1
                    sum_drivergene_hotspot+=1
                    keygene_dict[item_name].append(start)
                    
                if gene in CCG_list:
                   item[lencol+6]=item[lencol+6]+gene+', '
                   item[lencol+7]=item[lencol+7]+1
                   sum_cancercensusgene_hotspot+=1
                   keygene_dict[item_name].append(start)
    
    #print(gene_dict)               
    genetitle=['gene','num.gene','all cancer driver gene','num.all cancer driver gene',cancer+' driver gene','num.'+cancer+' driver gene','cancer census gene','num.cancer census gene']
    title=title+genetitle
    all.insert(0,title)    
    
    gene_statistic_table=['gene',sum_allGene_hotspot,sum_allGene_genome-sum_allGene_hotspot]
    statistic_table.loc[len(statistic_table)]=gene_statistic_table
    
    gene_statistic_table=['cancer census gene',sum_cancercensusgene_hotspot, sum_cancercensusgene_genome-sum_cancercensusgene_hotspot]
    statistic_table.loc[len(statistic_table)]=gene_statistic_table
    
    gene_statistic_table=[cancer+' driver gene',sum_drivergene_hotspot, sum_drivergene_genome-sum_drivergene_hotspot]
    statistic_table.loc[len(statistic_table)]=gene_statistic_table
    
    gene_statistic_table=['all cancer driver gene',sum_alldrivergene_hotspot,sum_alldrivergene_genome-sum_alldrivergene_hotspot]
    statistic_table.loc[len(statistic_table)]=gene_statistic_table
    
    
    

    ratio_allGene_genome=(sum_allGene_genome-sum_allGene_hotspot)/(3000000000-hotspotLenth)*1000000
    ratio_alldrivergene_genome=(sum_alldrivergene_genome-sum_alldrivergene_hotspot)/(3000000000-hotspotLenth)*1000000
    ratio_drivergene_genome=(sum_drivergene_genome-sum_drivergene_hotspot)/(3000000000-hotspotLenth)*1000000
    ratio_cancercensusgene_genome=(sum_cancercensusgene_genome-sum_cancercensusgene_hotspot)/(3000000000-hotspotLenth)*1000000
    
    ratio_allGene_hotspot=sum_allGene_hotspot/hotspotLenth*1000000
    ratio_alldrivergene_hotspot=sum_alldrivergene_hotspot/hotspotLenth*1000000
    ratio_drivergene_hotspot=sum_drivergene_hotspot/hotspotLenth*1000000
    ratio_cancercensusgene_hotspot=sum_cancercensusgene_hotspot/hotspotLenth*1000000
    
    #print('热点区的基因密度'+str(ratio_allGene_hotspot)+'\MB')
    #print('背景区的基因密度'+str(ratio_allGene_genome)+'\MB')
    #print('热点区的全部癌症驱动基因密度'+str(ratio_alldrivergene_hotspot)+'\MB')
    #print('背景区的全部癌症驱动基因密度'+str( ratio_alldrivergene_genome)+'\MB')
    #print('热点区的指定癌症驱动基因密度'+str(ratio_drivergene_hotspot)+'\MB')
    #print('背景区的指定癌症驱动基因密度'+str(ratio_drivergene_genome)+'\MB')
    #print('热点区的cancer census基因密度'+str(ratio_cancercensusgene_hotspot)+'\MB')
    #print('背景区的cancer census基因密度'+str(ratio_cancercensusgene_genome)+'\MB')
    
    return cancer+' driver gene',gene_dict,keygene_dict 
    #return ratio_allGene_hotspot,ratio_allGene_genome,ratio_alldrivergene_hotspot,ratio_alldrivergene_genome,ratio_drivergene_hotspot,ratio_drivergene_genome,ratio_cancercensusgene_hotspot,ratio_cancercensusgene_genome,gene_dict,keygene_dict