# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 15:54:29 2020

@author: Lv
"""

#Rscript code/findhotspot/findhotspot.R -i Results/signature3/SV.csv -k 8 -f signature3 -c 0

import os
import sys

sys.path.append('code/findhotspot')
import plothot




import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--method", help="hotspot detect method, first/seconde ,the first is fixed bins, the second is pcf", default='second')  
parser.add_argument("-f","--find", help="type/category/signature need to find hotspot,seprated by  ,") 
parser.add_argument("-k", "--kmin", help="smallest sample in every hotspot interval", default='8')  
parser.add_argument("-c", "--cluster", help="0:nonclustered,1:clustered", default='')  

args = parser.parse_args()

method=args.method
find=args.find
kmin=args.kmin
cluster=args.cluster

if method=='second':
    findlist=find.split(',')
    for item in findlist:
        command='Rscript code/findhotspot/findhotspot.R -i Results/'+item+'/SV.csv'+' -k '+str(kmin)+' -f '+item
        if cluster=='0' or cluster=='1':
            command=command+' -c '+cluster
        os.system(command)
        plothot.plothot(item)
        
            

