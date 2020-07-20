# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 10:18:03 2019

@author: Lv
"""

#https://www.jianshu.com/p/00425f6c0936  argparse
import os

#os.chdir('C:/Users/Lv/Desktop/annotation')
#sys.path.append('code/annotation')

import csv
import setSuperEnhancer
import setGene
import gwasCatalog
import eQTL
import plotfigure
import enhancer
import matplotlib.pyplot as plt
import pandas as pd 
import dataprocess
import annotother

def annotation_one_type(cancer,typefilename):
    print(typefilename+' annotating..')

    path='Results/'+typefilename+'/'
    
    start_list=[]
    
    with open(path+typefilename+'_hotspot.csv','r') as f:
    #with open('Results/signature3/hotspot_signature3.csv','r') as f:
        reader = csv.reader(f)
        for i in reader:
            start_list.append(i)
            
    title=start_list[0]
    del start_list[0]#先删掉列名一会再加上
    hotspotLenth=0
    
    for item in start_list:
        item[1]=int(item[1])
        item[2]=int(item[2])
        item[3]=int(item[3])
        hotspotLenth+=item[3]
    start_list.insert(0,title)
    
    
    
    
    #注释的列表
    #print(start_list[0])
    hotspot_id_index=len(start_list[0])-1
    
    
    
    statistic_table=pd.DataFrame(columns=['element name','num in hotspot','num in nonhotspot'])
    #statistic_table=pd.DataFrame(columns=['elementname','special type','other type','special type in hotspot','special type in nonhotspot','all type in hotspot','all type in nonhotspot'])
    
    
    
    #superenhancer_dict是为hotspot_match_sv服务的
    #ratio_allSEA_genome,ratio_allSEA_hotspot,ratio_SEA_genome,ratio_SEA_hotspot,superenhancer_dict= setSuperEnhancer.setSEA(start_list,hotspotLenth,hotspot_id_index)
    superenhancertitle,superenhancer_dict=setSuperEnhancer.setSEA(start_list, hotspot_id_index, statistic_table, cancer)
    #print(statistic_table)
    drivergenetitle,gene_dict,keygene_dict=setGene.setGene(start_list, hotspot_id_index, statistic_table, cancer,hotspotLenth)
    gwastitle,gwas_dict=gwasCatalog.setGWAS(start_list,hotspot_id_index, statistic_table, cancer,hotspotLenth)
    enhancertitle,enhancer_dict=enhancer.setEnhancer(start_list,hotspot_id_index, statistic_table, cancer,hotspotLenth)
    eQTLtitle=eQTL.seteQTL(start_list,hotspot_id_index, statistic_table, cancer,hotspotLenth)
    
    
    dataprocess.statistic_data_process(statistic_table,path,hotspotLenth)#处理表格
    
    #dataprocess.printlen(statistic_table,path)
    
    
    file=open(path+'annotation_'+typefilename+'.csv','w')#,newline='')
    testfile=csv.writer(file)
    
    #把注释的列表和本来的拼接在一起
    
    for i in range(len(start_list)):
    
        testfile.writerow(start_list[i])
    
    file.close()
    print('done') 
    
    plotfigure.plot_interval_percent(path,path+'annotation_'+typefilename+'.csv',superenhancertitle,enhancertitle,drivergenetitle,eQTLtitle,gwastitle)
    plotfigure.plot_distribution(path,path+'annotation_'+typefilename+'.csv',superenhancertitle,enhancertitle,gwastitle, eQTLtitle,drivergenetitle)
    plotfigure.plotratio(statistic_table,path)#画ratio图
    
    annotother.annotother(typefilename,hotspotLenth)
    
    
    
    

'''
cancer='LUAD'#胃癌
typefilename='lung'

annotation_one_type(cancer,typefilename)
'''

#给path和文件名
#(path,path+'annotation_'+typefilename+'.csv',superenhancertitle,enhancertitle,drivergenetitle,eQTLtitle,gwastitle)


'''
ratio_allGene_hotspot,ratio_allGene_genome,ratio_alldrivergene_hotspot,ratio_alldrivergene_genome,ratio_drivergene_hotspot,ratio_drivergene_genome,ratio_cancercensusgene_hotspot,ratio_cancercensusgene_genome\
,gene_dict,keygene_dict=setGene.setGene(start_list,hotspotLenth,hotspot_id_index)

ratio_GWAS_hotspot,ratio_GWAS_genome=gwasCatalog.setGWAS(start_list,hotspotLenth)
ratio_eQTL_hotspot,ratio_eQTL_genome=eQTL.seteQTL(start_list,hotspotLenth)
ratio_enhancer_hotspot,ratio_enhancer_genome=enhancer.setSEA(start_list,hotspotLenth)

#file=open('Results/signature3/annotation_signature3_2(enhancer).csv','w',newline='')
file=open('Results/lung/lung_annotation_2.csv','w',newline='')
#file=open('Results/liver/liver_annotation_2.csv','w',newline='')

#file=open('Results/annotation_gastric_nonclustered_del_8.csv','w',newline='')
#file=open('Results/rs3_1_annotation.csv','w',newline='')
#file=open('Results/prostate_annotation_27.csv','w',newline='')

testfile=csv.writer(file)

#把注释的列表和本来的拼接在一起

for i in range(len(start_list)):

    testfile.writerow(start_list[i])

file.close()
print('done') 

'''
#plotfigure.plotfigures(ratio_allSEA_genome,ratio_allSEA_hotspot,ratio_SEA_genome,ratio_SEA_hotspot,ratio_allGene_hotspot,ratio_allGene_genome,ratio_alldrivergene_hotspot,ratio_alldrivergene_genome,ratio_drivergene_hotspot,ratio_drivergene_genome,ratio_cancercensusgene_hotspot,ratio_cancercensusgene_genome,ratio_GWAS_hotspot,ratio_GWAS_genome,ratio_eQTL_hotspot,ratio_eQTL_genome)




