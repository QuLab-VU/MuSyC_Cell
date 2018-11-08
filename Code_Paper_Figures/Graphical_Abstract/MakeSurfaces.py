#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 14:03:53 2018

@author: xnmeyer
"""

import pandas as pd
import numpy as np
import matplotlib.pylab as plt
#from matplotlib.patches import FancyArrowPatch
from matplotlib import rc
rc('text', usetex=False)
font = {'family' : 'arial',
        'weight':'normal',
        'size'   : 8}
axes = {'linewidth': 2}
rc('font', **font)
rc('axes',**axes)
from plotting_surfaces import plot_surface
from plotting_surfaces import get_params
from plotting_surfaces import dip
from pythontools.synergy import synergy_tools, synergy_plots
from scipy.interpolate import interp1d
from scipy.optimize import root


fit_params = {}
fit_params['min_conc_d1'] = np.power(10.,-10)
fit_params['min_conc_d2'] = np.power(10.,-10)
fit_params['max_conc_d1'] = np.power(10.,-5)
fit_params['max_conc_d2'] = np.power(10.,-5)
fit_params['h1'] = 1.
fit_params['h2'] = 1.
fit_params['log_C1'] = -7.5
fit_params['log_C2'] = -6.5
fit_params['E0'] = .3
fit_params['E1'] = -0.1
fit_params['E2'] = 0.15
fit_params['E3'] = -0.1
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
zero_conc=2
zlim=(-0.35,.35)
plot_surface(fit_params,metric_name="Drug Effect",N=N,zero_conc=zero_conc,zlim=zlim,elev=20, azim=24,fname='original')



fit_params['E3'] = -0.3
fit_params['log_alpha1']=0.
fit_params['log_alpha2']=0.
plot_surface(fit_params,metric_name="Drug Effect",N=N,zero_conc=zero_conc,zlim=zlim,elev=20, azim=24,fname='Syn_eff')

fit_params['E3'] = -0.1
fit_params['log_alpha1']=2.
fit_params['log_alpha2']=2.
plot_surface(fit_params,metric_name="Drug Effect",N=N,zero_conc=zero_conc,zlim=zlim,elev=20, azim=24,fname='Syn_pot')
