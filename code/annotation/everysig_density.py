# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 10:55:20 2020
把之前注释的main.py 改成一个函数，解决一次性注释多个的问题
然后生成一个文件，是每个signature的密度合集，列名是signature1，signature2.。行名是SEA,allSEA...
@author: Lv
"""
import os
import sys
import pandas as pd
os.chdir('C:/Users/Lv/Desktop/annotation')
sys.path.append('code/annotation')

import os
import csv
import setSuperEnhancer
import setGene
import gwasCatalog
import eQTL
import plotfigure
import matplotlib.pyplot as plt

#name:如signature1,del_clustered_10kb,del
def annotation(name):
    start_list=[]
    #with open('Results/signature3/hotspot_signature3.csv','r') as f:
    with open('Results/'+name+'/hotspot_'+name+'.csv','r') as f:
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
    print(start_list[0])
    hotspot_id_index=len(start_list[0])-1
    
    #superenhancer_dict是为hotspot_match_sv服务的
    ratio_allSEA_genome,ratio_allSEA_hotspot,ratio_SEA_genome,ratio_SEA_hotspot,superenhancer_dict= setSuperEnhancer.setSEA(start_list,hotspotLenth,hotspot_id_index)
    
    ratio_allGene_hotspot,ratio_allGene_genome,ratio_alldrivergene_hotspot,ratio_alldrivergene_genome,ratio_drivergene_hotspot,ratio_drivergene_genome,ratio_cancercensusgene_hotspot,ratio_cancercensusgene_genome\
    ,gene_dict,keygene_dict=setGene.setGene(start_list,hotspotLenth,hotspot_id_index)
    
    ratio_GWAS_hotspot,ratio_GWAS_genome=gwasCatalog.setGWAS(start_list,hotspotLenth)
    ratio_eQTL_hotspot,ratio_eQTL_genome=eQTL.seteQTL(start_list,hotspotLenth)
    
    #if(not os.path.exists('Results/8')):
    #    os.makedirs('Results/8') 
    file=open('Results/'+name+'/annotation_'+name+'.csv','w',newline='')
    
    testfile=csv.writer(file)
    
    
    
    for i in range(len(start_list)):
    
        testfile.writerow(start_list[i])
    
    file.close()
    print('done') 
    
    plotfigure.plotratio(name,ratio_allSEA_genome,ratio_allSEA_hotspot,ratio_SEA_genome,ratio_SEA_hotspot,ratio_allGene_hotspot,ratio_allGene_genome,ratio_alldrivergene_hotspot,ratio_alldrivergene_genome,ratio_drivergene_hotspot,ratio_drivergene_genome,ratio_cancercensusgene_hotspot,ratio_cancercensusgene_genome,ratio_GWAS_hotspot,ratio_GWAS_genome,ratio_eQTL_hotspot,ratio_eQTL_genome)
    return ratio_allSEA_genome,ratio_allSEA_hotspot,ratio_SEA_genome,ratio_SEA_hotspot,ratio_allGene_hotspot,ratio_allGene_genome,ratio_alldrivergene_hotspot,ratio_alldrivergene_genome,ratio_drivergene_hotspot,ratio_drivergene_genome,ratio_cancercensusgene_hotspot,ratio_cancercensusgene_genome,ratio_GWAS_hotspot,ratio_GWAS_genome,ratio_eQTL_hotspot,ratio_eQTL_genome


ratiodata=pd.DataFrame(index=['all super enhancer','related super enhancer','gene','cancer census gene','all cancer driver gene','related driver gene','related GWAS SNP','related eQTL'])
for i in range(9):
   ratio_allSEA_genome,ratio_allSEA_hotspot,ratio_SEA_genome,ratio_SEA_hotspot,ratio_allGene_hotspot,ratio_allGene_genome,ratio_alldrivergene_hotspot,ratio_alldrivergene_genome,ratio_drivergene_hotspot,ratio_drivergene_genome,ratio_cancercensusgene_hotspot,ratio_cancercensusgene_genome,ratio_GWAS_hotspot,ratio_GWAS_genome,ratio_eQTL_hotspot,ratio_eQTL_genome=annotation('signature'+str(i+1))
   ratiodata['signature'+str(i+1)]=[round(ratio_allSEA_hotspot/ratio_allSEA_genome,2),round(ratio_SEA_hotspot/ratio_SEA_genome,2),round(ratio_allGene_hotspot/ratio_alldrivergene_genome,2),round(ratio_cancercensusgene_hotspot/ratio_cancercensusgene_genome,2),round(ratio_alldrivergene_hotspot/ratio_alldrivergene_genome,2),round(ratio_drivergene_hotspot/ratio_drivergene_genome,2),round(ratio_GWAS_hotspot/ratio_GWAS_genome,2),round(ratio_eQTL_hotspot/ratio_eQTL_genome,2)]
ratiodata.to_csv('Results/ratio.csv')
#data = pd.DataFrame(index=['a','b'],columns=['s','d'])
#data['c'] = ['d','2']
#data.to_csv("tzzs_data.csv")