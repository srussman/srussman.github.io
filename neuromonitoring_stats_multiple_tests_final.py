#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 10:24:28 2021

@author: srussman@ucsd.edu

This code runs statistical tests from the paper
"""

## Load data and packages
import pandas as pd
import os
import scipy.stats as stats
import numpy as np
import statsmodels.stats.multitest as smt

#os.chdir('Downloads/')

###### Multiple tests ######

# %% Two-sample, paired t-test
# Question: how to account for difference in variance?
# Assumptions: equal variance, normal distribution
df = pd.read_excel("subdural1_LU_30_snr.xlsx") # import file
df.columns = [c.replace(' ', '') for c in df.columns]
df=df.transpose()
array = pd.DataFrame(df).to_numpy()
array=np.concatenate( array, axis=0 )

df2 = pd.read_excel("subdural1_LU_30_clean.xlsx")
N=int(len(df2)/2)
indices1 = np.arange(0,N)
indices2 = np.arange(N,2*N)
idx1=indices1[array]
idx2=indices2[array]


res = stats.ttest_rel(df2['d_1'][idx1], df2['d_1'][idx2])
display(res)

# %% Mann Whitney U test (Wilcoxon paired test)

res2 = stats.wilcoxon(df2['d_1'][idx1], df2['d_1'][idx2])
display(res2)

# %% Kruskal Wallis test

res3 = stats.kruskal(df2['d_1'][idx1], df2['d_1'][idx2])
display(res3)

# %% Calculate Cohen's d

M1 = np.mean(df2['d_1'][idx1])
M2 = np.mean(df2['d_1'][idx2])
SD1 = np.std(df2['d_1'][idx1])
SD2 = np.std(df2['d_1'][idx2])

sd_pooled = np.sqrt((SD1**2 + SD2**2)/2)

Cohens_d = (M2 - M1)/sd_pooled * (N-3)/(N-2.25)*np.sqrt((N-2)/N)
Cohens_d_corr = (M2 - M1)/sd_pooled
display(Cohens_d)
display(Cohens_d_corr)


# Save data
res_df = pd.DataFrame(res)
res2_df = pd.DataFrame(res2)
res3_df = pd.DataFrame(res3)
concat_df = pd.concat([res_df,res2_df,res3_df])
cohen_df = pd.DataFrame({'Cohen': [Cohens_d]})
cohen_corr_df = pd.DataFrame({'Cohen': [Cohens_d_corr]})
cohen_concat_df = pd.concat([cohen_df,cohen_corr_df])

concat_df.to_csv("subdural1_LU_30_tests_sample.csv")
cohen_concat_df.to_csv("subdural1_LU_30_cohen_sample.csv")

#%% Run Whtmey Mann U test on each channel and save
import os
import glob
dat_f = []
filelist=[]
mydir = 'Downloads/'
for path, subdirs, files in os.walk(mydir):
    for file in files:
        if (file.endswith('manova_data.xlsx') & file.startswith('subdural1_LU_75')):
            filelist.append(os.path.join(path, file))
from tkinter import Tcl
filelist=Tcl().call('lsort', '-dict', filelist)

for i in range(len(filelist)): #LU$: 367, LU15: 367,733, LU30: 733,1098, LU75: 1098,1464
    file=filelist[i]
    display(file)
    df2=pd.read_excel(file)
                  
    full_diff = []
    for ch in range(len(df2.columns)):
        dat_seg = df2[df2.columns[ch]]
        mi=min(dat_seg)
        ma=max(dat_seg)
        diff = ma-mi
        full_diff.append(diff)

    df = pd.DataFrame(full_diff)
    idx = int(len(df2.columns)/2)
    res = stats.wilcoxon(df[0][0:idx], df[0][idx:(idx*2)])

    res_df = pd.DataFrame(res)
    dat_f.append(res_df[0][1])
    display(res_df)

dat_f_df = pd.DataFrame(dat_f)
dat_f_df.to_csv("subdural1_LU_75_mannwhit.csv")


