#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 14:03:53 2018

@author: xnmeyer
"""

import pandas as pd
import numpy as np
from plotting_surfaces import plot_surface
import matplotlib.pyplot as plt
import matplotlib.cm as cm

fit_params = {}
fit_params['min_conc_d1'] = np.power(10.,-10)
fit_params['min_conc_d2'] = np.power(10.,-10)
fit_params['max_conc_d1'] = np.power(10.,-5)
fit_params['max_conc_d2'] = np.power(10.,-5)
fit_params['h1'] = 1.
fit_params['h2'] = 1.
fit_params['log_C1'] = -7.5
fit_params['log_C2'] = -6.5
fit_params['E0'] = 1
fit_params['E1'] = .5
fit_params['E2'] = .5
fit_params['E3'] = 0.
fit_params['log_alpha1']=0.
fit_params['log_alpha2']=0.
fit_params['log_gamma1']=0.
fit_params['log_gamma2']=0.
fit_params['r1']=100000.*fit_params['E0']
fit_params['r2']=100000.*fit_params['E0']
fit_params['drug1_name']='d1'
fit_params['drug2_name']='d2'
fit_params['expt'] = 'cartoon'
fit_params = pd.DataFrame([fit_params])
N=100
zero_conc=1
zlim=(0,1)


fig = plt.figure(figsize=(10,10),facecolor='k')
ax=[];e3=np.linspace(0,.35,12)[::-1];a=np.linspace(-1.5,1.5,13);
al=[];e=[]
for i in range(13):
    if np.remainder(i,2)==0:
        for j in range(11):
            ax.append(fig.add_axes([0.06+.08*(j), 0.9-.07*(i), 0.1, 0.1],projection='3d'))        
            al.append(a[i])
            e.append(e3[j+1])
    else:
        for j in range(12):
            ax.append(fig.add_axes([0.02+.08*(j), 0.9-.07*(i), 0.1, 0.1],projection='3d'))
            al.append(a[i])
            e.append(e3[j])

for i in range(len(ax)):  
    fit_params['E3'] = e[i]
    fit_params['log_alpha1']=al[i]
    fit_params['log_alpha2']=al[i]
    plot_surface(fit_params,ax=ax[i],cmap=cm.viridis,incl_gamma=False, N=50, elev=25, azim=45, zlim=zlim, zero_conc=1, alpha = 1)



plt.savefig('cover_illusion.png',bbox_inches = 'tight',pad_inches = 0)



