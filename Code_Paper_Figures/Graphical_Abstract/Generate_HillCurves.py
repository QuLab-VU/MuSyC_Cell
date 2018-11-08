#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 11:26:52 2018

@author: xnmeyer
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
font = {'family' : 'normal',
        'weight':'normal',
        'size'   : 8}
axes = {'linewidth': 1}
rc('font', **font)
rc('axes',**axes)


def Edrug1D(d,h,E0,Em,C):
    return Em + (E0-Em) / (1 + (d/C)**h)


        
####For figure 1
fig = plt.figure(figsize = (1,2.5))
d1 = np.logspace(-10,-5,1000);
h1 = 1;
C1 = 1e-7;
E0 = 1;
E1 = 0;

Ed = np.zeros((1000,5))
ax1 = plt.subplot(312)
for e,h in enumerate([.5,1,2,10]):
    Ed[:,e] = Edrug1D(d1,h,E0,E1,C1)    
    plt.plot(np.log10(d1),Ed[:,e],color='k',linewidth=1)
ax1.spines['right'].set_color('none')
ax1.spines['top'].set_color('none')
ax1.spines['bottom'].set_color('none')
plt.axhline(0, color='k')
plt.xticks([])
plt.yticks([])   


Ed = np.zeros((1000,5))
ax2 = plt.subplot(311)
for e,Em in enumerate([.45,.3,.15,0]):
    Ed[:,e] = Edrug1D(d1,h1,E0,Em,C1)    
    plt.plot(np.log10(d1),Ed[:,e],color='k',linewidth=1)
ax2.spines['right'].set_color('none')
ax2.spines['top'].set_color('none')
ax2.spines['bottom'].set_color('none')
plt.axhline(0, color='k')
plt.xticks([])
plt.yticks([])  


Ed = np.zeros((1000,5))
ax3 = plt.subplot(313)
for e,C in enumerate([1e-6,1e-7,1e-8,1e-9]):
    Ed[:,e] = Edrug1D(d1,h,E0,E1,C)    
    plt.plot(np.log10(d1),Ed[:,e],color='k',linewidth=1)
ax3.spines['right'].set_color('none')
ax3.spines['top'].set_color('none')
ax3.spines['bottom'].set_color('none')
plt.axhline(0, color='k')
plt.xticks([])
plt.yticks([])  
plt.tight_layout()


plt.savefig('hill_Curves.pdf')
