#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:11:38 2017

@author: xnmeyer
"""
import pandas as pd
import numpy as np
from itertools import compress
import matplotlib.pyplot as plt
from matplotlib import rc
font = {'family' : 'normal',
        'weight':'normal',
        'size'   : 8}
axes = {'linewidth': 2}
rc('font', **font)
rc('axes',**axes)
from scipy.integrate import quad
from matplotlib.artist import setp
import itertools
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec


##############################################################
##Read in the data and filter out the NSCLC cells to leave only the melanoma cells
###############################################################
T = pd.read_csv('../../Data/MasterResults.csv')
T = T[T['cell_line']!='PC9C1']
T = T[T['cell_line']!='BR1']
T = T[T['cell_line']!='PC9AZR']
T = T[T['MCMC_converge']==1]
T = T[T['model_level']>4]
T = T[T['R2']>.5]

###############################################################
#Make waterfall plots and scatter plots
###############################################################
c = plt.get_cmap('nipy_spectral')
cols = c(np.linspace(0,1,9))
a = []
for letter in range(97,97+16):
    a.append(chr(letter))
m = []
for i in range(16):
    m.append('$' + a[i] + '$')
    
    
fig = plt.figure(figsize=(7.5,5))
ax = []
for i in range(16):
    if i<4:
        v = 0 
    elif i<8:
        v = 1 
    elif i<12:
        v = 2 
    else:
        v = 3
    if i<4:
        l = i 
    elif i<8:
        l = i-4 
    elif i<12:
        l = i-8 
    else:
        l = i-12
    ax.append(plt.subplot2grid((4,5),(v,l)))
ax.append(plt.subplot2grid((4,5),(0,4),rowspan=3))
plt.subplots_adjust(wspace=0,hspace=.4)


cell_lines = list(T.drop_duplicates('cell_line')['cell_line'])
t = T.drop_duplicates(['drug1_name','drug2_name'])[['drug1_name','drug2_name']]
drugs = list(np.sort(t['drug1_name'] + '_' + t['drug2_name']))
p = []
seen_d = [];seen_c = [];idx_c = [];idx_d = [];patch = [];mark = [];keep=[];
for e,i in enumerate(T.index):
    idx1 = cell_lines.index(T['cell_line'].loc[i])
    idx2 = drugs.index(T['drug1_name'].loc[i] + '_' + T['drug2_name'].loc[i])
    plt.sca(ax[idx1+8])
    p.append(plt.scatter(T['log_alpha1'].loc[i],T['log_alpha2'].loc[i],s=30,c=cols[idx1,:],marker=m[idx2],label=T['drug1_name'].loc[i] + '_' + T['drug2_name'].loc[i]))
    plt.title(T['cell_line'].loc[i])
    if seen_d.count(idx2)==0:
        seen_d.append(idx2)
        patch.append(mpatches.Patch(color=cols[idx1,:], label=T['drug1_name'].loc[i] + '_' + T['drug2_name'].loc[i]))
        idx_d.append(T['drug1_name'].loc[i] + '_' + T['drug2_name'].loc[i])
        keep.append(e)
    if seen_c.count(idx1)==0:
        seen_c.append(idx1)
        mark.append(T['cell_line'].loc[i])
        idx_c.append(e)
        
x_min = np.min([np.min(T[['log_alpha1','log_alpha2']],axis=1)])-2
x_max = np.max([np.max(T[['log_alpha1','log_alpha2']],axis=1)])+2
y_min = np.min([np.min(T[['log_alpha1','log_alpha2']],axis=1)])-2
y_max = np.max([np.max(T[['log_alpha1','log_alpha2']],axis=1)])+2
for e in range(8):
    ax[e+8].set_xlim((x_min,x_max))
    ax[e+8].set_ylim((y_min,y_max))
    y_lim = (y_min,y_max)
    x_lim = (x_min,x_max)
    plt.sca(ax[e+8])
    plt.plot([0,0],y_lim,linestyle='--',color='k')
    plt.plot(x_lim,[0,0],linestyle='--',color='k')
    ax[e+8].tick_params(direction='in')
    if e+8!=12:
        ax[e+8].set_xticklabels([])
        ax[e+8].set_yticklabels([])
        ax[e+8].set_xticks([])
plt.sca(ax[12])    
plt.xlabel(r'$log(\alpha_1)$')
plt.ylabel(r'$log(\alpha_2)$')
plt.sca(ax[-1])
ax[12].tick_params(direction='out')
pp = [p[i] for i in keep]
ppt = pp[:]
for e,i in enumerate(seen_d):
    ppt[e]=pp[np.where(np.in1d(seen_d,e))[0][0]]
legend1 = plt.legend(handles=ppt,loc=2)
for i in range(len(legend1.legendHandles)):
    legend1.legendHandles[i].set_color('k')
ax[-1].axis('off')

y_max = max(T['beta_obs_norm'])+.5
y_min = min(T['beta_obs_norm'])-.5

T.reset_index(drop = True,inplace=True)
T['drug_label'] = 'hi'
for i in T.index:
    T.set_value(i,'drug_label',T['drug1_name'].loc[i] + '_' + T['drug2_name'].loc[i])

for e in range(8):
    plt.sca(ax[e])
    df_tmp = T[T['cell_line'] ==  mark[e]]
    idx = []
    for i in df_tmp['drug_label']:
        idx.append(drugs.index(i))
    tmp_m = [m[i] for i in idx] 
    y = np.array(df_tmp['beta_obs_norm'])
    y_sd=np.array(df_tmp['beta_obs_norm_std'])
    x = np.argsort(y)
    plt.bar(range(len(x)),[y[i] for i in x],color=cols[seen_c[e],:])
    plt.errorbar(range(len(x)),[y[i] for i in x],yerr=[y_sd[i] for i in x],color='k',fmt=' ')

#    ax[e].set_xticks(range(len(x)))
#    ax[e].set_xticklabels([tmp_m[i] for i in x],va='top')
    tmp_m = [tmp_m[i] for i in x]
    y = [y[i] for i in x]
    x = range(len(x))
    for i in x:
        if y[i]>0:
            plt.annotate(tmp_m[i],xy = (i-.5,y[i]+.1))
        else:
            plt.annotate(tmp_m[i],xy = (i-.5,y[i]-.8))
    ax[e].set_xticks([])
    if e!=0 and e!=4:
        ax[e].tick_params(direction='in')
        ax[e].set_yticklabels([])
    else:
        plt.ylabel(r'$\beta_{obs}$')
    ax[e].set_ylim((y_min-.5,y_max+.5))
    plt.title(mark[e])

plt.savefig('Plots/Figure3_DSDs.pdf')


