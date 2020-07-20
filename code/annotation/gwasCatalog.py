# -*- coding: utf-8 -*-
import indextable
def setGWAS(start_list,hotspot_id_index, statistic_table, cancer,hotspotLenth):
    
    type=indextable.indextable(cancer ,'gwas')
    
    all=start_list
    title=all[0]
    del all[0]
    
    lencol=len(all[0])
    for item in all: 
       
        item.append(0)
        item.append('')
        item.append(0)
        
    
    sum_GWAS_genome=0
    sum_GWAS_hotspot=0 
    sum_allGWAS_genome=0   
    sum_allGWAS_hotspot=0    
    
    gwas_dict={}        
    for item in all:
        item_name=item[hotspot_id_index]
        gwas_dict[item_name]=[]
    
    #f=open('data/gwas.txt','r',encoding='UTF-8')
    f=open('data/gwas.txt','r')
 
    for line in f:
        sum_allGWAS_genome+=1
        if type in line or type+' cancer' in line:

            sum_GWAS_genome+=1
      
        linelist=line.split('\t')
        #disease=line[7]
        chr=linelist[11].split(';')[0]#为什么是分号，因为有的就是3；4；5；
        chr=chr.split(' ')[0]
        pos=linelist[12].split(';')[0]#有的就是12765097；3443452
        pos=pos.split(' ')[0]
        if chr=='':
            continue
        if pos=='':
            continue 
        snp=linelist[21]
        posnum=int(pos)
  
        for item in all:
            item_name=item[hotspot_id_index]
            if (chr==item[0]) and (posnum>item[1] and posnum<item[2]):
                sum_allGWAS_hotspot+=1
                item[lencol]=item[lencol]+1
                if type in line or type+' cancer' in line:
                    
                    item[lencol+1]=item[lencol+1]+snp+', '
                    item[lencol+2]=item[lencol+2]+1
                    sum_GWAS_hotspot+=1
                    gwas_dict[item_name].append(posnum)
                  
    snptitle=['num.all gwas snp',type+' cancer snp','num.'+type+' cancer snp']
    title=title+snptitle
    all.insert(0,title)

    ratio_GWAS_genome=(sum_GWAS_genome-sum_GWAS_hotspot)/(3000000000-hotspotLenth)*1000000
    ratio_GWAS_hotspot=sum_GWAS_hotspot/hotspotLenth*1000000   
    #print('指定癌症的snp在热点区的密度是'+ str(ratio_GWAS_hotspot)+'/MB')
    #print('指定癌症的snp在背景区的密度是'+ str(ratio_GWAS_genome)+'/MB')    
    
    cellbase_list=['all gwas snp',sum_allGWAS_hotspot,sum_allGWAS_genome-sum_allGWAS_hotspot]
    statistic_table.loc[len(statistic_table)]=cellbase_list
    
    cellbase_list=[type+' cancer related gwas snp', sum_GWAS_hotspot, sum_GWAS_genome-sum_GWAS_hotspot]
    statistic_table.loc[len(statistic_table)]=cellbase_list
    
    return type+' cancer snp',gwas_dict