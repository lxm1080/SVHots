# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 19:03:34 2019

@author: Lv
"""



import indextable


def seteQTL(start_list,hotspot_id_index, statistic_table, cancer,hotspotLenth):
    eQTLfilelist,type=indextable.indextable(cancer,'eQTL')
    if(len(eQTLfilelist)==0):
        return
    
    all=start_list
    title=all[0]
    del all[0]
    
    lencol=len(all[0])
    for item in all: 
        item.append('')
        item.append(0)
        item.append('')
        item.append(0)
    
    sum_eQTL_genome=0
    sum_eQTL_hotspot=0    
    
    for i in range(len(eQTLfilelist)): 
        file=open('data/eQTL/'+eQTLfilelist[i]+'.signifpairs.txt','r')
       
        for line in file.readlines()[1:len(file.readlines())-1]:
            sum_eQTL_genome+=1
            line=line.split('\t')
            variant_id=line[0]
            egene_id=line[1]
            
            chr=variant_id.split('_')[0][3:]
            pos=variant_id.split('_')[1]
            pos=int(pos)
        
            for item in all:
                if (chr==item[0]) and (pos>item[1] and pos<item[2]):
         #           item[lencol]=item[lencol]+variant_id+', '
                    item[lencol+1]=item[lencol+1]+1
                    sum_eQTL_hotspot+=1    
                    if egene_id not in item[lencol+2]:
                        item[lencol+2]=item[lencol+2]+egene_id+', '
                        item[lencol+3]=item[lencol+3]+1
                


    
    eQTLtitle=[type+' related eQTL','num.'+type+' related eQTL',type + ' related egene','num.'+type +' related egene']
    title=title+eQTLtitle
    all.insert(0,title)
    
    ratio_eQTL_genome=(sum_eQTL_genome-sum_eQTL_hotspot)/(3000000000-hotspotLenth)*1000000
    ratio_eQTL_hotspot=sum_eQTL_hotspot/hotspotLenth*1000000
    
    #print('热点区的指定组织eQTL的密度是'+str(ratio_eQTL_hotspot)+'/MB')
    #print('背景区的指定组织eQTL的密度是'+str(ratio_eQTL_genome)+'/MB')
    
    #cellbase_list=['eQTL',type+' related eQTL', 'other eQTL', sum_eQTL_hotspot, sum_eQTL_genome-sum_eQTL_hotspot,5000000,6000000]
    cellbase_list=[type+' related eQTL', sum_eQTL_hotspot, sum_eQTL_genome-sum_eQTL_hotspot]
    statistic_table.loc[len(statistic_table)]=cellbase_list
    #return ratio_eQTL_hotspot,ratio_eQTL_genome
    return type+' related eQTL'