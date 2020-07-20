# -*- coding: utf-8 -*-
"""
Created on Fri May 22 17:40:08 2020

@author: Lv
"""



import os
import sys
import annotationOneType
import plotfigure

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-a", "--annot", help="type/signatue/category name to annotate,seprated by ',' ", default='',action="store_true")   # action 为"store_true"或"store_false"
parser.add_argument("--cancer", help="cancer type of samples", default='STAD',action="store_true")   # action 为"store_true"或"store_false"
parser.add_argument("-c", "--compare", help="type/signatue/category name to compare,seprated by ',' ", default='',action="store_true")   # action 为"store_true"或"store_false"
args = parser.parse_args()

cancer=args.cancer
annotationStr=args.a
compareStr=args.c


#os.chdir('C:/Users/Lv/Desktop/annotation')
sys.path.append('code/annotation')

#输入
#cancer='STAD'
#annotationStr='signature1,signature2,signature3,signature4,signature5,signature6,signature7,signature8,signature9'
#compareStr='signature1,signature2,signature3,signature4,signature5,signature6,signature7,signature8,signature9'






if annotationStr != '':
    annotationlist=annotationStr.split(',')
    for i in range(len(annotationlist)):    
        typefilename=annotationlist[i]
        annotationOneType.annotation_one_type(cancer,typefilename)
        

if compareStr != '':
    comparelist=compareStr.split(',')
    plotfigure.plot_OR_comparison(comparelist)




#typelist=['signature1','signature2','signature3','signature4','signature5','signature6','signature7','signature8','signature9',]
