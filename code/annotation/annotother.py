# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 19:58:58 2020
annot others
@author: Lv
"""


import pandas as pd
import os
import sys
sys.path.append('code/annotation')
import dataprocess


def annotcount(typefilename,annotname,annotation):
    #typefilename='signature3'
    #annotname='dgv'
    print('startcount')
    annotfile='data/'+annotname+'.bed'
    hotspot_bed='Results/'+typefilename+'/hotspot.bed'
    command='bedtools intersect -a '+hotspot_bed +' -b '+annotfile+  ' -wa -wb | bedtools groupby -i - -g 1-4 -c 8 -o collapse >'+'Results/'+typefilename+'/'+annotname+'_results.bed'
    os.system(command)
    #和dgv的bed合并
    dgv_dict={}
    dgv_annot=open('Results/'+typefilename+'/'+annotname+'_results.bed','r')
    for line in dgv_annot:
        line=line.split('\t')
        dgv_dict[line[3]]=line[4]
    dgv_annot.close()
    annotation[annotname]=['']*len(annotation.index)
    for key in dgv_dict.keys():
        annotation.loc[[key],[annotname]]=len(dgv_dict[key].split(',')) #注释个数
        #annotation.loc[[key],[annotname]]=dgv_dict[key]#注释列表
    print('endcount')


def annotfiles(annotation,filestype,typefilename,hotspot_len):
    print('startfiles')
    #filestype='TAD'
    files = os.listdir('data/'+filestype)
    hotspot_bed='Results/'+typefilename+'/hotspot.bed'
    statistic_table=pd.DataFrame(columns=['element name','num of intervals','num in hotspot','num in nonhotspot','hotspot length','nonhotspot length','p(one sided poisson test)'])
    
    for file in files:#对每一个bed文件    
        annotfile='data/'+filestype+'/'+file 
        numingenome=int(os.popen('cat '+  annotfile+'|wc -l').read())
        command='bedtools intersect -b '+hotspot_bed +' -a '+annotfile+  ' -wa -wb | bedtools groupby -i - -g 1-4 -c 8 -o collapse | wc -l'
        numinhot=os.popen(command)
        numinhot=int(numinhot.read())
        print(file+':'+str(numinhot))
       
        
        elementname=file.split('.')[0]
        numinnonhot=numingenome-numinhot
        len_nonhot=3000-hotspot_len
        p=dataprocess.my_poisson(numinhot,hotspot_len,numinnonhot,len_nonhot)
        
        #算一下一共多少热点区间出现了这个element
        command1='bedtools intersect -a '+hotspot_bed +' -b '+annotfile+  ' -wa -wb | bedtools groupby -i - -g 1-4 -c 8 -o collapse | wc -l'
        interval_num=os.popen(command1)
        interval_num=interval_num.read()    
        
        newlist=[elementname, interval_num , numinhot, numinnonhot, hotspot_len, len_nonhot,p]
        statistic_table.loc[len(statistic_table)]=newlist
    print('endfiles')
    
    
    return statistic_table  





def annotother(typefilename,hotspotLenth):
    
    #hotspot_len=100000000
    hotspot_len=round(hotspotLenth/1000000,3)
    
    #typefilename='signature3'
    #根据typefilename找到annotation
    annotation=pd.read_csv('Results/'+typefilename+'/'+'annotation_'+typefilename+'.csv')
    chr=annotation['chr']
    start=annotation['start.bp']
    end=annotation['end.bp']
    id=annotation['hotspot.id']
    annotation.index=annotation['hotspot.id']
    #写hotspot.bed文件
    hotspot_bed=open('Results/'+typefilename+'/hotspot.bed','w')
    for i in range(len(chr)):
        hotspot_bed.write('chr'+str(chr[i]))
        hotspot_bed.write('\t')
        hotspot_bed.write(str(start[i]))
        hotspot_bed.write('\t')
        hotspot_bed.write(str(end[i]))
        hotspot_bed.write('\t')
        hotspot_bed.write(str(id[i]))
        #hotspot_bed.write('\t')
        hotspot_bed.write('\n')
        i+=1
    hotspot_bed.close()
    
    annotcount(typefilename,'dgv',annotation)
    
    
    TADtable=annotfiles(annotation,'TAD',typefilename,hotspot_len)
    TFBStable=annotfiles(annotation,'TFBS',typefilename,hotspot_len)
    
    
    
    #def fileline():
     
    #annotation.to_csv('Results/'+typefilename+'/'+'annotation_'+typefilename+'.csv',index=False) #保存，并且不保存索引项
    with pd.ExcelWriter('Results/'+typefilename+'/'+'annotation_'+typefilename+'.xlsx') as writer:
          annotation.to_excel(writer, sheet_name='annotation',index=False)
          TADtable.to_excel(writer, sheet_name='TAD',index=False)
          TFBStable.to_excel(writer, sheet_name='TFBS',index=False)
