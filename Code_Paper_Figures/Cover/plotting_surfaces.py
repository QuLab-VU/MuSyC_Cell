#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 08:59:05 2018

Code to plot the dose response surface in python matplotlib
Data points are optional
Auto-scaled to surface limits
Colormap is scaled to surface unless specified
Based on code from David Wooten
"""


####################################################################
###Import Packages
####################################################################
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import *
from matplotlib.collections import PolyCollection
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import colorConverter
#from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
from matplotlib import rc
from scipy import interpolate
rc('text', usetex=False)
font = {'family' : 'arial',
        'weight':'normal',
        'size'   : 8}
axes = {'linewidth': 3}
rc('font', **font)
rc('axes',**axes)


##########################################
### D E F I N E   F U N C T I O N S
##########################################    
#incl_gamma=True; data_pts=None; N=50; elev=20; azim=19; fname=None; zlim=(-.1,1.1); zero_conc=1; alpha = 0.9;metric_name="per eff";path=path = 'conflating_plot/Malarial_combinations/';title=None

def plot_surface(fit_params,ax=None,incl_gamma=False, N=50, elev=20, azim=19, zlim=None,cmap=cm.viridis, zero_conc=1, alpha = 0.9):

    e0, e1, e2, e3, h1, h2, r1, r1r, r2, r2r, ec50_1, ec50_2, alpha1, alpha2, gamma1,gamma2,d1_min, d1_max, d2_min, d2_max, drug1_name, drug2_name, expt = get_params(fit_params)
    if ax is None:
        fig = plt.figure(figsize=(3.5,3))
        ax = fig.gca(projection='3d')
    #ax.set_axis_off()
    y = np.linspace(d1_min-zero_conc, d1_max ,N)
    t = np.linspace(d2_min-zero_conc, d2_max ,N)

    if incl_gamma:
        tt, yy = np.meshgrid(t, y)
        conc3 = np.array([(np.power(10.,c1),np.power(10.,c2)) for c1 in y for c2 in t])
        zz = np.array([Edrug2D_NDB_hill(d,e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2) for d in conc3])
        zz = zz.reshape(N,N)
    else:
        yy, tt = np.meshgrid(y, t)
        zz = dip(yy,tt, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
    if zlim is None:
        zmin = np.min(zz)
        zmax = np.max(zz)
        if np.abs(zmin) > np.abs(zmax): zmax = np.abs(zmin)
        else: zmin = -np.abs(zmax)
    else:
        zmin = zlim[0]
        zmax = zlim[1]
        
    # Plot the surface for real, with appropriate alpha
    surf = ax.plot_surface(tt, yy, zz, cstride=1, rstride=1, alpha=alpha, cmap = cmap, vmin=zmin, vmax=zmax, linewidth=0)
    
    # colored curves on left and right
    lw = 1
    btop_map = np.ones((N,4))
    btop_map[:,1]=0
    btop_map[:,0]=np.linspace(0,1,N) 
    rtop_map = np.ones((N,4))
    rtop_map[:,1]=0
    rtop_map[:,2]=np.linspace(0,1,N)    
    if incl_gamma:
        ax.plot(d2_min*np.ones(y.shape)-zero_conc, y, Edrug2D_NDB_hill((np.power(10.,y),np.power(10.,d2_min*np.ones(y.shape)-zero_conc)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), linewidth=lw,color='k')
        ax.scatter(d2_max*np.ones(y.shape), y, Edrug2D_NDB_hill((np.power(10.,y),np.power(10.,d2_max*np.ones(y.shape))),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), s=lw*10,c='k')

        ax.plot(t, d1_min*np.ones(y.shape)-zero_conc, Edrug2D_NDB_hill((np.power(10.,d1_min*np.ones(t.shape)-zero_conc),np.power(10.,t)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), linewidth=lw,color='k')
        ax.scatter(t, d1_max*np.ones(y.shape), Edrug2D_NDB_hill((np.power(10.,d1_max*np.ones(t.shape)),np.power(10.,t)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), s=lw*10,c='k')
    
    else:
        ax.plot(d2_min*np.ones(y.shape)-zero_conc, y, dip(y,d2_min-zero_conc, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), linewidth=lw,color='k')
        ax.scatter(d2_max*np.ones(y.shape), y, dip(y,d2_max, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), s=lw/10.,c='k')
    
        ax.plot(t, d1_min*np.ones(y.shape)-zero_conc, dip(d1_min-zero_conc,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), linewidth=lw,color='k')
        ax.scatter(t, d1_max*np.ones(y.shape), dip(d1_max,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), s=lw/10.,c='k')
        
    
    if incl_gamma:
        # light grey grid across surface
        for ttt in np.linspace(d2_min-zero_conc,d2_max,10):
            ax.plot(ttt*np.ones(y.shape), y, Edrug2D_NDB_hill((np.power(10.,y),np.power(10.,ttt)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), '-k', linewidth=1, alpha=0.1)
            ax.plot(ttt*np.ones(y.shape), y, Edrug2D_NDB_hill((np.power(10.,y),np.power(10.,ttt)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), '-k', linewidth=1, alpha=0.1)
            ax.plot(ttt*np.ones(y.shape), y, Edrug2D_NDB_hill((np.power(10.,y),np.power(10.,ttt)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), '-k', linewidth=1, alpha=0.1)
    
        for yyy in np.linspace(d1_min-zero_conc, d1_max,10):
            ax.plot(t, yyy*np.ones(y.shape), Edrug2D_NDB_hill((np.power(10.,yyy),np.power(10.,t)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), '-k', linewidth=1, alpha=0.1)
            ax.plot(t, yyy*np.ones(y.shape), Edrug2D_NDB_hill((np.power(10.,yyy),np.power(10.,t)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), '-k', linewidth=1, alpha=0.1)
            ax.plot(t, yyy*np.ones(y.shape), Edrug2D_NDB_hill((np.power(10.,yyy),np.power(10.,t)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), '-k', linewidth=1, alpha=0.1)
        
    else:
        # light grey grid across surface
        for ttt in np.linspace(d2_min-zero_conc,d2_max,10):
            ax.plot(ttt*np.ones(y.shape), y, dip(y,ttt, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)
            ax.plot(ttt*np.ones(y.shape), y, dip(y,ttt, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)
            ax.plot(ttt*np.ones(y.shape), y, dip(y,ttt, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)
    
        for yyy in np.linspace(d1_min-zero_conc, d1_max,10):
            ax.plot(t, yyy*np.ones(y.shape), dip(yyy,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)
            ax.plot(t, yyy*np.ones(y.shape), dip(yyy,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)
            ax.plot(t, yyy*np.ones(y.shape), dip(yyy,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)

    # Set the view
    ax.view_init(elev=elev, azim=azim)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_aspect('equal')
    ax.patch.set_visible(False)

    

def get_params(df):
    #Effects
    e0 = np.float(df['E0'])
    e1 = np.float(df['E1'])
    e2 = np.float(df['E2'])
    e3 = np.float(df['E3'])
    #EC50s
    ec50_1 = np.power(10.,np.float(df['log_C1']))
    ec50_2 = np.power(10.,np.float(df['log_C2']))
    #Hill slopes
    h1 = np.float(df['h1'])
    h2 = np.float(df['h2'])
    #Transition parameters
    r1 = np.float(df['r1'])
    r2 = np.float(df['r2'])
    r1r = r1*np.power(ec50_1,h1)
    r2r = r2*np.power(ec50_2,h2)
    #synergy parameters
    alpha_1 = np.power(10.,np.float(df['log_alpha1']))
    alpha_2 = np.power(10.,np.float(df['log_alpha2']))
    
    #synergy parameters
    gamma_1 = np.power(10.,np.float(df['log_gamma1']))
    gamma_2 = np.power(10.,np.float(df['log_gamma2']))
    
    #max and min concentrations for the two drugs    
    concentration_1_min = np.log10(np.float(df['min_conc_d1']))
    concentration_2_min = np.log10(np.float(df['min_conc_d2']))
    concentration_1_max = np.log10(np.float(df['max_conc_d1']))
    concentration_2_max = np.log10(np.float(df['max_conc_d2']))
    #Get drug names
    drug1_name = str(df['drug1_name'].values[0])
    drug2_name = str(df['drug2_name'].values[0])
    #Get experiment name
    expt = str(df['expt'].values[0])
    return e0, e1, e2, e3, h1, h2, r1, r1r, r2, r2r, ec50_1, ec50_2, alpha_1, alpha_2, gamma_1, gamma_2, concentration_1_min, concentration_1_max, concentration_2_min, concentration_2_max, drug1_name, drug2_name, expt

#2D hill equation (full model...ie no detail balance)
def dip(d1, d2, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r):
    d1 = np.power(10.,d1);      d2 = np.power(10.,d2)
    h1 = 1.*h1;                 h2 = 1.*h2
    alpha1 = 1.*alpha1;         alpha2 = 1.*alpha2
    r1 = 1.*r1;                 r1r = 1.*r1r
    r2 = 1.*r2;                 r2r = 1.*r2r
    #Percent unaffected in each population
    U = r1r*r2r*(r1*(alpha2*d1)**h1 + r1r + r2*(alpha1*d2)**h2 + r2r)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)
    A1 = r1*r2r*(d1**h1*r1*(alpha2*d1)**h1 + d1**h1*r1r + d1**h1*r2r + d2**h2*r2*(alpha2*d1)**h1)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)
    A2 = r1r*r2*(d1**h1*r1*(alpha1*d2)**h2 + d2**h2*r1r + d2**h2*r2*(alpha1*d2)**h2 + d2**h2*r2r)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)
    A12 = r1*r2*(d1**h1*r1*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r2r*(alpha1*d2)**h2 + d2**h2*r1r*(alpha2*d1)**h1 + d2**h2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)
    #Return dip rate
    return U*e0 + A1*e1 + A2*e2 + A12*e3




def Edrug2D_NDB_hill(d,E0,E1,E2,E3,r1,r2,C1,C2,h1,h2,alpha1,alpha2,gamma1,gamma2):
    d1 = d[0]
    d2 = d[1]
    Ed =    (E2*alpha1**(gamma1*h2)*d2**(h2*(gamma1 + 1))*r2**(gamma1 + 1)*(C1**h1*r1)**gamma2 + 
             E1*alpha2**(gamma2*h1)*d1**(h1*(gamma2 + 1))*r1**(gamma2 + 1)*(C2**h2*r2)**gamma1 + 
             C1**h1*E0*r1**(gamma2 + 1)*(alpha2*d1)**(gamma2*h1)*(C2**h2*r2)**gamma1 + 
             C2**h2*E0*r2**(gamma1 + 1)*(alpha1*d2)**(gamma1*h2)*(C1**h1*r1)**gamma2 + 
             E2*d1**h1*r1**(gamma2 + 1)*r2**gamma1*(C1**h1)**gamma2*(alpha1*d2)**(gamma1*h2) + 
             E1*d2**h2*r1**gamma2*r2**(gamma1 + 1)*(C2**h2)**gamma1*(alpha2*d1)**(gamma2*h1) + 
             E3*alpha1**(gamma1*h2)*d2**(h2*(gamma1 + 1))*r1**gamma2*r2**(gamma1 + 1)*(alpha2*d1)**(gamma2*h1) + 
             E3*alpha2**(gamma2*h1)*d1**(h1*(gamma2 + 1))*r1**(gamma2 + 1)*r2**gamma1*(alpha1*d2)**(gamma1*h2) + 
             C1**h1*C2**h2*E0*r1*r2**(gamma1 + 1)*(C2**h2)**gamma1 + C1**h1*C2**h2*E0*r1**(gamma2 + 1)*r2*(C1**h1)**gamma2 + 
             C2**h2*E3*d1**h1*r1*r2**(gamma1 + 1)*(alpha1*d2)**(gamma1*h2) + 
             C1**h1*E3*d2**h2*r1**(gamma2 + 1)*r2*(alpha2*d1)**(gamma2*h1) + 
             C2**h2*E1*d1**h1*r1*r2**(gamma1 + 1)*(C2**h2)**gamma1 + 
             C2**h2*E1*d1**h1*r1**(gamma2 + 1)*r2*(C1**h1)**gamma2 + 
             C1**h1*E2*d2**h2*r1*r2**(gamma1 + 1)*(C2**h2)**gamma1 + 
             C1**h1*E2*d2**h2*r1**(gamma2 + 1)*r2*(C1**h1)**gamma2)/ \
            (C1**h1*r1**(gamma2 + 1)*(alpha2*d1)**(gamma2*h1)*(C2**h2*r2)**gamma1 + 
             C2**h2*r2**(gamma1 + 1)*(alpha1*d2)**(gamma1*h2)*(C1**h1*r1)**gamma2 + 
             alpha1**(gamma1*h2)*d2**(h2*(gamma1 + 1))*r2**(gamma1 + 1)*(C1**h1*r1)**gamma2 + 
             alpha2**(gamma2*h1)*d1**(h1*(gamma2 + 1))*r1**(gamma2 + 1)*(C2**h2*r2)**gamma1 + 
             alpha1**(gamma1*h2)*d2**(h2*(gamma1 + 1))*r1**gamma2*r2**(gamma1 + 1)*(alpha2*d1)**(gamma2*h1) + 
             alpha2**(gamma2*h1)*d1**(h1*(gamma2 + 1))*r1**(gamma2 + 1)*r2**gamma1*(alpha1*d2)**(gamma1*h2) + 
             C1**h1*C2**h2*r1*r2**(gamma1 + 1)*(C2**h2)**gamma1 + C1**h1*C2**h2*r1**(gamma2 + 1)*r2*(C1**h1)**gamma2 + 
             C2**h2*d1**h1*r1*r2**(gamma1 + 1)*(alpha1*d2)**(gamma1*h2) + 
             C1**h1*d2**h2*r1**(gamma2 + 1)*r2*(alpha2*d1)**(gamma2*h1) +
             C1**h1*d2**h2*r1*r2**(gamma1 + 1)*(C2**h2)**gamma1 +
             C1**h1*d2**h2*r1**(gamma2 + 1)*r2*(C1**h1)**gamma2 + 
             C2**h2*d1**h1*r1*r2**(gamma1 + 1)*(C2**h2)**gamma1 + 
             C2**h2*d1**h1*r1**(gamma2 + 1)*r2*(C1**h1)**gamma2 + 
             d1**h1*r1**(gamma2 + 1)*r2**gamma1*(C1**h1)**gamma2*(alpha1*d2)**(gamma1*h2) + 
             d2**h2*r1**gamma2*r2**(gamma1 + 1)*(C2**h2)**gamma1*(alpha2*d1)**(gamma2*h1))
    return Ed




