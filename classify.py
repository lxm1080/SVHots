# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 10:05:14 2019

@author: Lv
"""





import os
import sys
sys.path.append('code/classify')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



#画饼图
def func(pct, allvals):
    absolute = int(pct/100.*np.sum(allvals))
    return "{:.1f}%".format(pct, absolute)

def plotpie(method,typelist):
    #method='type'
    #typelist=['del','type','dup','inv','inv']
    
    type_dict={}
    for thistype in typelist:
        if  thistype in type_dict.keys():
            type_dict[thistype]+=1
        else:
            type_dict[thistype]=1
    
    typename=[]
    for name in type_dict.keys():
        typename.append(name)
        
    data=[]
    for val in type_dict.values():
        data.append(val)
       
    fig, ax = plt.subplots(figsize=(10, 5), subplot_kw=dict(aspect="equal"))
    
    ingredients = typename
    #wedges, texts, autotexts = ax.pie(data, autopct=lambda pct: func(pct, data),
    #                                  textprops=dict(color="w",fontsize=20))    
    wedges, texts, autotexts = ax.pie(data, autopct=lambda pct: func(pct, data),
                                      textprops={'color' : 'white', 'size': 20} )    
    leg=ax.legend(wedges, ingredients,
              #title=method,
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1),
              fontsize=15)
    leg.get_frame().set_linewidth(0.0)
    
    plt.setp(autotexts, size=8, weight="bold")
    
    #ax.set_title("distriubution of structural variants")
    plt.savefig('Results/classify.jpg',dpi=200, bbox_inches='tight')








import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--classify", help="classify method, no/type/category/signature", default='signature')  
parser.add_argument("-i","--input", help="input file, SV information") 

args = parser.parse_args()

method=args.classify
inputfile=args.input

#method='type'
#inputfile='liver_data.csv'


if os.path.exists('Results') :
    print()
else:
    os.mkdir("Results") 


SVfile= pd.read_csv(inputfile) #读入csv

if method == 'no':
    newfile='Results/no'
    os.mkdir(newfile)
    SVfile.to_csv(newfile+'/SV.csv',index=False) 
    
if method =='type':
    typelist=SVfile['Type']
    plotpie(method,typelist)
    typelist = list(set(typelist)) 
    for onetype in typelist:
        onetypeSV=SVfile.loc[SVfile["Type"] == onetype] #筛选得到Type是某个种类的
        newfile='Results/'+onetype
        os.mkdir(newfile)
        onetypeSV.to_csv(newfile+'/SV.csv',index=False) 

#会生成signature.csv到Results，生成classify文件夹到
if method =='category' or method =='signature':
    if os.path.exists('Results/signature.csv') :
        SVfile= pd.read_csv('Results/signature.csv')
    else:
        command='Rscript code/classify/testsignature.R --input '+inputfile
        os.system(command)
        SVfile= pd.read_csv('Results/signature.csv')

    
if method =='category':
    typelist=SVfile['Category']
    plotpie(method,typelist)
    typelist = list(set(typelist)) 
    for onetype in typelist:
        onetypeSV=SVfile.loc[SVfile["Category"] == onetype] #删选Type是某个值的
        newfile='Results/'+onetype
        os.mkdir(newfile)
        onetypeSV.to_csv(newfile+'/SV.csv',index=False) 
        
if method =='signature':
    SVfile=SVfile.sort_values(by="Sig.max" , ascending=True)
    typelist=SVfile['Sig.max']
    plotpie(method,typelist)
    typelist = list(set(typelist)) 
    for onetype in typelist:
        onetypeSV=SVfile.loc[SVfile[onetype+'.prob']>0.5] #删选Type是某个值的
        newfile='Results/'+onetype.replace(".","").lower()
        os.mkdir(newfile)
        onetypeSV.to_csv(newfile+'/SV.csv',index=False) 
     



    


