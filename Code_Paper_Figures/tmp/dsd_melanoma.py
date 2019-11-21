#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 10:54:51 2019

@author: meyerct6
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

import pandas as pd
import numpy as np
import sys
import os
from itertools import compress
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import rc
#Adjust text from  https://github.com/Phlya/adjustText
#Can be installed with pip
from adjustText import adjust_text
import scipy.stats as st
font = {'family' : 'normal',
        'weight':'normal',
        'size'   : 8}
axes = {'linewidth': 2}
rc('font', **font)
rc('axes',**axes)
from scipy import integrate
from scipy.interpolate import interp2d

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

num_classes = len(np.unique(T['cell_line']))
cell_lines = np.unique(T['cell_line'])

xmin,xmax = (T['log_alpha2'].min()-0.5,T['log_alpha2'].max()+0.5)
ymin,ymax = (T['beta_obs_norm'].min()-.05,T['beta_obs_norm'].max()+.1)
y_lim = (ymin,ymax)
x_lim = (xmin,xmax)

#Colors for each super class
col = cm.nipy_spectral(np.linspace(0,1,len(cell_lines)))

fig = plt.figure(figsize=(8.5,2.3))
ax = []
for i in range(num_classes):
    ax.append(plt.subplot2grid((1,num_classes),(0,i)))
plt.subplots_adjust(wspace=0,hspace=.4)

for i in range(num_classes):
    plt.sca(ax[i])
    x = T.loc[T['cell_line']==cell_lines[i],'log_alpha2'].values
    y = T.loc[T['cell_line']==cell_lines[i],'beta_obs_norm'].values
    plt.scatter(x,y,s=30,color=col[i],edgecolor='k',linewidth=1)
    plt.title(cell_lines[i])
    
for e in range(num_classes):
    ax[e].set_xlim((xmin,xmax))
    ax[e].set_ylim((ymin,ymax))
    plt.sca(ax[e])
    plt.plot([0,0],y_lim,linestyle='--',color='k')
    plt.plot(x_lim,[0,0],linestyle='--',color='k')
    plt.sca(ax[e])
    #plt.legend(loc=9, bbox_to_anchor=(0.5, -0.03), ncol=1)
    ax[e].tick_params(direction='in')
    ax[e].set_xticklabels([])
    ax[e].set_yticklabels([])
    ax[e].xaxis.set_ticks_position('bottom')
    ax[e].yaxis.set_ticks_position('left')
    
plt.sca(ax[0])
plt.xlabel(r'$log(\alpha_2)$'+'\n(MEKi potentiates BRAFi)',labelpad=-4)
plt.ylabel(r'$\bar{\beta_{obs}}$')
