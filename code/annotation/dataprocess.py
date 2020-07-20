# -*- coding: utf-8 -*-
"""
Created on Sun May 10 13:56:29 2020
处理数据
@author: Lv
"""
import os
#import sys
#os.chdir('C:/Users/Lv/Desktop/annotation')
#sys.path.append('code/annotation')
import pandas as pd




#from decimal import *
#getcontext().prec = 6
#print(Decimal(1)/Decimal(7))
'''
def my_kf(a,b,c,d):
    kf_data = np.array([[a,b], [c,d]])
    kf = chi2_contingency(kf_data)
    pvalue=kf[1]
    getcontext().prec = 3
    return Decimal(pvalue)/Decimal(1)
'''
def my_poisson(a,b,c,d):
    #a= speicial type in hot, b= hot len, c= special type in nonhot , d=nonhotlen
    command='Rscript code/annotation/poisson_test.R'+' --a '+str(a)+' --b '+str(b)+' --c '+str(c)+' --d '+str(d)
    
    p=os.popen(command) #在linux中运行时换成这个
    p=p.read()
    #print('p:'+p)
    
    #p='2.4e-19'
    #getcontext().prec = 3
    return float(p)
    
def my_minus(a, b):
    return a-b
def my_divide(a,b):
    if b==0:
        return 0
    else:
        return round(a/b,4)



#这里开始def处理df的函数
def statistic_data_process(df,path,hotspot_len):
    
    hotspot_len=round(hotspot_len/1000000,3)
    df['hotspot region length']=[hotspot_len]*len(df)
    df['nonhotspot region length']=[3000-hotspot_len]*len(df)
    
    df['density in hotspot']=df.apply(lambda col: my_divide(col['num in hotspot'],col['hotspot region length']), axis=1)
    df['density in nonhotspot']=df.apply(lambda col: my_divide(col['num in nonhotspot'],col['nonhotspot region length']), axis=1)
    df['OR']=df.apply(lambda col: my_divide(col['density in hotspot'],col['density in nonhotspot']), axis=1)
    df['pvalue(poisson test)']=df.apply(lambda col: my_poisson(col['num in hotspot'],col['hotspot region length'],col['num in nonhotspot'], col['nonhotspot region length']), axis=1)

       
    df.to_csv(path+'cancer_specifity_statistic.csv', sep=',', header=True, index=True)
    

'''
nonhotspot_len=2850
hotspot_len=150

df=pd.DataFrame(columns=['elementname','special type','other type','special type in hotspot','special type in nonhotspot','all type in hotspot','all type in nonhotspot'])

#每个函数里加一加
row1=['super_enhancer','breast_related super_enhancer','other super enhancer',2,10,5,100]
row2=['gene','cancer_census_gene','other gene',2,20,3,100]

#处理一下这个表格
df.loc[len(df)]=row1
df.loc[len(df)]=row2
'''


'''
#这里开始def处理df的函数
def statistic_data_process(df,path,hotspot_len):
    print(len(df))
    hotspot_len=hotspot_len/1000000
    df['other type in hotspot'] = df.apply(lambda col: my_minus(col['all type in hotspot'], col['special type in hotspot']), axis=1)
    df['other type in nonhotspot'] = df.apply(lambda col: my_minus(col['all type in nonhotspot'], col['special type in nonhotspot']), axis=1)
    df['hotspot region length']=[hotspot_len]*len(df)
    df['nonhotspot region length']=[3000-hotspot_len]*len(df)

    df['special type density in hotspot']=df.apply(lambda col: my_divide(col['special type in hotspot'],col['hotspot region length']), axis=1)
    df['special type density in nonhotspot']=df.apply(lambda col: my_divide(col['special type in nonhotspot'],col['nonhotspot region length']), axis=1)
    df['OR of special type']=df.apply(lambda col: my_divide(col['special type density in hotspot'],col['special type density in nonhotspot']), axis=1)
    df['pvalue_special type']=df.apply(lambda col: my_poisson(col['special type in hotspot'],col['hotspot region length'],col['special type in nonhotspot'], col['nonhotspot region length']), axis=1)
    
    df['other type density in hotspot']=df.apply(lambda col: my_divide(col['other type in hotspot'],col['hotspot region length']), axis=1)
    df['other type density in nonhotspot']=df.apply(lambda col: my_divide(col['other type in hotspot'],col['nonhotspot region length']), axis=1)
    df['OR of other type']=df.apply(lambda col: my_divide(col['other type density in hotspot'],col['other type density in nonhotspot']), axis=1)
    df['pvalue_other type']=df.apply(lambda col: my_poisson(col['other type in hotspot'],col['hotspot region length'],col['other type in nonhotspot'], col['nonhotspot region length']), axis=1)

    
    df['all type density in hotspot']=df.apply(lambda col: my_divide(col['all type in hotspot'],col['hotspot region length']), axis=1)
    df['all type density in nonhotspot']=df.apply(lambda col: my_divide(col['all type in nonhotspot'],col['nonhotspot region length']), axis=1)
    df['OR of all type']=df.apply(lambda col: my_divide(col['all type density in hotspot'],col['all type density in nonhotspot']), axis=1)
    df['pvalue_all type']=df.apply(lambda col: my_poisson(col['all type in hotspot'],col['hotspot region length'],col['all type in nonhotspot'], col['nonhotspot region length']), axis=1)

    
    df['OR of special type and other type ']=df.apply(lambda col: my_divide(col['OR of special type'],col['OR of other type']), axis=1)
    
    df.to_csv(path+'cancer_specifity_statistic.csv', sep=',', header=True, index=True)

#statistic_data_process(df,'Results/signatre3')函数没问题啊
#传入df画图，以下都是画图的代码
'''


