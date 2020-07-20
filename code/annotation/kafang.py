# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 19:46:45 2020

@author: Lv
"""

import numpy as np
import pandas as pd
from  scipy.stats import chi2_contingency
kf_data = np.array([[15,95], [85,5]])
kf = chi2_contingency(kf_data)
pvalue=kf[1]/2
#print(pvalue)
#print(kf)
#print('chisq-statistic=%.4f, p-value=%.4f, df=%i expected_frep=%s'%kf)


from scipy import stats

import numpy as np
#https://blog.csdn.net/u010891397/article/details/83823822
A = np.array([ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
B = np.array([ 0, 2, 5, 6, 8, 10, 13, 14, 17, 20])

stats.ttest_ind(B,A,equal_var= False)
print(stats.ttest_ind(B,A,equal_var= False))