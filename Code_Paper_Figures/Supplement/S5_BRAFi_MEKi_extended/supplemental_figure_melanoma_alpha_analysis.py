#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 11:25:11 2017

@author: xnmeyer
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
font = {'family' : 'normal',
        'weight':'normal',
        'size'   : 8}
axes = {'linewidth': 2}
rc('font', **font)
rc('axes',**axes)
np.random.seed(12345) #set a random seed so the mcmc fits are reproducible


T = pd.read_csv('../../../Data/MasterResults.csv')
T = T[T['cell_line']!='PC9C1']
T = T[T['cell_line']!='BR1']
T = T[T['cell_line']!='PC9AZR']

###############################################################
#Make waterfall plots and scatter plots
###############################################################
T = T[T['MCMC_converge']==1]
T = T[T['model_level']>4]

c = plt.get_cmap('spectral')
cols = ['xkcd:red purple','tab:blue', 'tab:orange', 'tab:red', 'tab:purple','tab:green','xkcd:blue green','tab:brown']


#Set up figure dimensions
fit = plt.figure(figsize=(7.5,9))
ax = []
for i in range(4):
    ax.append(plt.subplot2grid((7,2),(i,0)))
for i in range(4):
    ax.append(plt.subplot2grid((7,2),(i,1)))
ax.append(plt.subplot2grid((7,2),(4,0),colspan=2,rowspan=2))
ax.append(plt.subplot2grid((7,2),(6,0),colspan=2))
plt.subplots_adjust(wspace=0.2,hspace=.40)


t = T.groupby(['drug1_name', 'drug2_name'])
for e,keys in enumerate(t.groups.keys()):
    ax_idx = np.where(np.unique(T['drug1_name'])==keys[0])[0][0]
    x_idx =  np.where(np.unique(T['drug2_name'])==keys[1])[0][0]
    plt.sca(ax[ax_idx])
    y = t['log_alpha1'].get_group(keys).astype('float')
    cl = t['cell_line'].get_group(keys)
    pt = plt.boxplot(np.array(y),vert=True,widths=.5,showfliers=False,positions=[x_idx])
    x = x_idx+np.random.uniform(-.1,.1,len(y))
    for i,c in enumerate(cl):
        v = np.where(c==np.unique(T['cell_line']))[0][0]
        plt.scatter(x[i],np.array(y)[i],s=20,color=cols[v])
    plt.title(keys[0])
    ax[ax_idx].set_xlim((-.5,3.5))
    ax[ax_idx].set_ylim((T['log_alpha1'].min()-1,T['log_alpha2'].max()+1))
    if ax_idx == 3:
        ax[ax_idx].set_xticks([0,1,2,3])
        ax[ax_idx].set_xticklabels(np.unique(T['drug2_name']))
    else:
        ax[ax_idx].set_xticks([0,1,2,3])
        ax[ax_idx].set_xticklabels([])
    if ax_idx != 0:
        ax[ax_idx].set_yticks([])
    else:
        plt.ylabel(r'$log(\alpha_1)$')
    plt.plot([-.5,3.5],[0,0],'k--')
    
for e,keys in enumerate(t.groups.keys()):
    ax_idx = np.where(np.unique(T['drug2_name'])==keys[1])[0][0]
    x_idx =  np.where(np.unique(T['drug1_name'])==keys[0])[0][0]
    plt.sca(ax[ax_idx+4])
    y = t['log_alpha2'].get_group(keys).astype('float')
    cl = t['cell_line'].get_group(keys)
    pt = plt.boxplot(np.array(y),vert=True,widths=.5,showfliers=False,positions=[x_idx])
    x = x_idx+np.random.uniform(-.1,.1,len(y))
    for i,c in enumerate(cl):
        v = np.where(c==np.unique(T['cell_line']))[0][0]
        plt.scatter(x[i],np.array(y)[i],s=20,c=cols[v])
    plt.title(keys[1])
    ax[ax_idx+4].set_xlim((-.5,3.5))
    ax[ax_idx+4].set_ylim((T['log_alpha1'].min()-1,T['log_alpha2'].max()+1))
    if ax_idx+4 == 7:
        ax[ax_idx+4].set_xticks([0,1,2,3])
        ax[ax_idx+4].set_xticklabels(np.unique(T['drug1_name']))
    else:
        ax[ax_idx+4].set_xticks([0,1,2,3])
        ax[ax_idx+4].set_xticklabels([])
    if ax_idx+4 != 4:
        ax[ax_idx+4].set_yticks([])
    else:
        plt.ylabel(r'$log(\alpha_2)$')
    plt.plot([-.5,3.5],[0,0],'k--')
    
    
###############################################################
#Make boxplots for beta
###############################################################
t = T.groupby(['drug1_name', 'drug2_name'])
key = 'beta_obs_norm'
m = []
for e,keys in enumerate(t.groups.keys()):
    m.append(np.median(t[key].get_group(keys).astype('float')))
idx = np.argsort(m)

plt.sca(ax[-2])
for e,keys in enumerate(t.groups.keys()):
    y = t[key].get_group(keys).astype('float')
    cl = t['cell_line'].get_group(keys)
    pt = plt.boxplot(np.array(y),vert=True,widths=.5,showfliers=False,positions=[np.where(e==idx)[0]])
    x = np.where(e==idx)[0]+np.random.uniform(-.01,.01,len(y))
    for i,c in enumerate(cl):
        jit = np.random.uniform(low=-.2,high=.2)
        v = np.where(c==np.unique(T['cell_line']))[0][0]
        plt.scatter(x[i]+jit,np.array(y)[i],s=20,c=cols[v])
plt.plot([-1,len(t.groups.keys())],[0,0],'k--',)
ax[-2].set_xlim((-.6,len(t.groups.keys())-.4))
plt.ylabel(r'$\beta_{obs}$')
ax[-2].set_xticks(np.linspace(0,len(t.groups.keys())-1,len(t.groups.keys())))
v = t.groups.keys()
v = [v[i] for i in idx]
v = [k[0]+'\n'+k[1] for k in v]
ax[-2].set_xticklabels(v)
plt.xticks(rotation=90)

    
plt.sca(ax[-1])
for i,e in enumerate(np.unique(T['cell_line'])):
    plt.annotate(e,xy=(i,-.2),color=cols[i],ha='center',va='center',fontsize=8,fontweight='bold',annotation_clip=False)
ax[-1].set_xlim((-.5,7.5))
ax[-1].set_ylim((-.2,.2))
ax[-1].axis('off')
    
plt.savefig('SuppMelanoma_Alpha_Beta_plot.pdf')    
    
