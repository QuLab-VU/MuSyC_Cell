#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Christian T Meyer
email christian.t.meyer@vanderbilt.edu
Code for making both static and dynamic surface plots
in paper 'Quantifying Drug Synergy with respect to potency and efficacy'
Commented and up to date 4-27-18
"""


####################################################################
###Use David's code to create the surfaces:
####################################################################
#Import packages
import pandas as pd
import numpy as np
from pylab import *
from mpl_toolkits.mplot3d import *
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
from scipy import interpolate
rc('text', usetex=False)
font = {'family' : 'normal',
        'weight':'normal',
        'size'   : 8}
axes = {'linewidth': 2}
rc('font', **font)
rc('axes',**axes)
import plotly.graph_objs as go
from plotly.offline import  plot
import os
##########################################
### D E F I N E   F U N C T I O N S
##########################################
def get_params(df):
    e0 = np.float(df['E0'])
    e1 = np.float(df['E1'])
    e2 = np.float(df['E2'])
    e3 = np.float(df['E3'])

    h1 = np.float(df['h1'])
    h2 = np.float(df['h2'])

    r1 = np.float(df['r1'])
    r2 = np.float(df['r2'])

    ec50_1 = np.power(10.,np.float(df['log_C1']))
    ec50_2 = np.power(10.,np.float(df['log_C2']))

    r1r = r1*np.power(ec50_1,h1)
    r2r = r2*np.power(ec50_2,h2)

    alpha_1 = np.power(10.,np.float(df['log_alpha1']))
    alpha_2 = np.power(10.,np.float(df['log_alpha2']))
    beta = np.float(df['beta'])/(e0-np.min((e1,e2)))
    beta_obs = np.float(df['beta_obs'])/(e0-df[['E1_obs','E2_obs']].min())

    
    concentration_1_min = np.log10(np.float(df['min_conc_d1']))
    concentration_2_min = np.log10(np.float(df['min_conc_d2']))
    
    concentration_1_max = np.log10(np.float(df['max_conc_d1']))
    concentration_2_max = np.log10(np.float(df['max_conc_d2']))
    
    model_level = str(df['model_level'])
    
    drug1_name = str(df['drug1_name'])
    drug2_name = str(df['drug2_name'])

    expt = str(df['expt'])
    
    return e0, e1, e2, e3, h1, h2, r1, r1r, r2, r2r, ec50_1, ec50_2, alpha_1, alpha_2, beta, beta_obs, concentration_1_min, concentration_1_max, concentration_2_min, concentration_2_max, model_level, drug1_name, drug2_name, expt

def dip(d1, d2, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r):
#    d1 = 1.*d1
#    d2 = 1.*d2
    d1 = np.power(10.,d1)
    d2 = np.power(10.,d2)
    h1 = 1.*h1
    h2 = 1.*h2
    alpha1 = 1.*alpha1
    alpha2 = 1.*alpha2
    r1 = 1.*r1
    r1r = 1.*r1r
    r2 = 1.*r2
    r2r = 1.*r2r
    
    U = r1r*r2r*(r1*(alpha2*d1)**h1 + r1r + r2*(alpha1*d2)**h2 + r2r)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)

    A1 = r1*r2r*(d1**h1*r1*(alpha2*d1)**h1 + d1**h1*r1r + d1**h1*r2r + d2**h2*r2*(alpha2*d1)**h1)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)

    A2 = r1r*r2*(d1**h1*r1*(alpha1*d2)**h2 + d2**h2*r1r + d2**h2*r2*(alpha1*d2)**h2 + d2**h2*r2r)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)

    A12 = r1*r2*(d1**h1*r1*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r2r*(alpha1*d2)**h2 + d2**h2*r1r*(alpha2*d1)**h1 + d2**h2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)
    
    return U*e0 + A1*e1 + A2*e2 + A12*e3
    
def plot_surface(mcmc_data, hts_data, drug_id, drug, N=50, elev=20, azim=19, fname=None, zlim=None,zero_conc=1):
    e0, e1, e2, e3, h1, h2, r1, r1r, r2, r2r, ec50_1, ec50_2, alpha1, alpha2, beta, beta_obs, d1_min, d1_max, d2_min, d2_max, model_level, drug1_name, drug2_name, expt = get_params(mcmc_data)
    fig = figure(figsize=(11,6))
    ax = fig.gca(projection='3d')
    #ax.set_axis_off()
    y = linspace(d1_min-zero_conc, d1_max ,N)
    t = linspace(d2_min-zero_conc, d2_max ,N)
    yy, tt = meshgrid(y, t)
    zz = dip(yy,tt, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
    zmin = np.min(zz)
    zmax = np.max(zz)

    if np.abs(zmin) > np.abs(zmax): zmax = np.abs(zmin)
    else: zmin = -np.abs(zmax)

    # Plot it once on a throwaway axis, to get cbar without alpha problems
    my_cmap_rgb = plt.get_cmap('bwr')(np.arange(256))
    alpha = 0.8

    for i in range(3): # Do not include the last column!
        my_cmap_rgb[:,i] = (1 - alpha) + alpha*my_cmap_rgb[:,i]
    my_cmap = ListedColormap(my_cmap_rgb, name='my_cmap')
    surf = ax.plot_surface(tt, yy, zz, cstride=1, rstride=1, cmap = my_cmap, vmin=zmin, vmax=zmax, linewidth=0)
    cbar = fig.colorbar(surf)
    cbar.solids.set_rasterized(True)
    cbar.solids.set_edgecolor("face")

    cla()
    # Plot the surface for real, with appropriate alpha
    surf = ax.plot_surface(tt, yy, zz, cstride=1, rstride=1, alpha=alpha, cmap = cm.bwr, vmin=zmin, vmax=zmax, linewidth=0)
    
    # colored curves on left and right
    lw = 5
    ax.plot(d2_min*ones(y.shape)-zero_conc, y, dip(y,d2_min-zero_conc, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), linewidth=lw)
    ax.plot(d2_max*ones(y.shape), y, dip(y,d2_max, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), linewidth=lw)

    ax.plot(t, d1_min*ones(y.shape)-zero_conc, dip(d1_min-zero_conc,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), linewidth=lw)
    ax.plot(t, d1_max*ones(y.shape), dip(d1_max,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), linewidth=lw)
    
    
    # light grey grid across surface
    for ttt in linspace(d2_min-zero_conc,d2_max,10):
        ax.plot(ttt*ones(y.shape), y, dip(y,ttt, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)
        ax.plot(ttt*ones(y.shape), y, dip(y,ttt, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)
        ax.plot(ttt*ones(y.shape), y, dip(y,ttt, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)

    for yyy in linspace(d1_min-zero_conc, d1_max,10):
        ax.plot(t, yyy*ones(y.shape), dip(yyy,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)
        ax.plot(t, yyy*ones(y.shape), dip(yyy,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)
        ax.plot(t, yyy*ones(y.shape), dip(yyy,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), '-k', linewidth=1, alpha=0.1)

    # Set the view
    ax.view_init(elev=elev, azim=azim)


    ax.set_ylabel("log(%s)[M]"%mcmc_data['drug1_name'].values[0])
    ax.set_xlabel("log(%s)[M]"%mcmc_data['drug2_name'].values[0])
    ax.set_zlabel("DIP rate" + r'$h^{-1}$')


    # Scatter points, and get error_bars
    scat_d1 = np.log10(hts_data['drug1.conc'])
    scat_d1.loc[scat_d1==-np.inf] = scat_d1.loc[scat_d1!=-np.inf].min()-zero_conc
    scat_d2 = np.log10(hts_data['drug2.conc'])
    scat_d2.loc[scat_d2==-np.inf] = scat_d2.loc[scat_d2!=-np.inf].min()-zero_conc
    scat_dip = hts_data['rate']
    scat_erb = hts_data['rate.95ci']
    ax.scatter(scat_d2, scat_d1, scat_dip, s=5, depthshade=True)

    # Plot error bars
    for _d1, _d2, _dip, _erb in zip(scat_d1, scat_d2, scat_dip, scat_erb):
        ax.plot([_d2,_d2], [_d1,_d1], [_dip-_erb, _dip+_erb], 'k-', alpha=0.3,linewidth=1)


    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])

    # Plane of DIP=0
    c_plane = colorConverter.to_rgba('k', alpha=0.2)
    verts = [array([(d2_min-zero_conc,d1_min-zero_conc), (d2_min-zero_conc,d1_max), (d2_max,d1_max), (d2_max,d1_min-zero_conc), (d2_min-zero_conc,d1_min-zero_conc)])]
    poly = PolyCollection(verts, facecolors=c_plane)
    ax.add_collection3d(poly, zs=[0], zdir='z')

    # Plot intersection of surface with DIP=0 plane
    CS = contour(tt,yy,zz,levels=[0], linewidths=3, colors='k')

    # If needed, manually set zlim
    if zlim is not None: ax.set_zlim(zlim)

    # Save plot, or show it
    if fname is None: show()
    else: plt.savefig(fname)
    clf()
    plt.close()

def plot_slices(mcmc_data, hts_data, drug_id, drug, N=50, fname=None, lw=1,zero_conc=1):
    e0, e1, e2, e3, h1, h2, r1, r1r, r2, r2r, ec50_1, ec50_2, alpha1, alpha2, beta, beta_obs, d1_min, d1_max, d2_min, d2_max, model_level, drug1_name, drug2_name, expt = get_params(mcmc_data)
    fig = figure(figsize=(1.5,1))
    y = linspace(d1_min-zero_conc, d1_max, N)
    t = linspace(d2_min-zero_conc, d2_max, N)    
    dip1_0 = dip(y,-10**zero_conc, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
    dip1_1 = dip(y,d2_max, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
    dip2_0 = dip(-10**zero_conc,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
    dip2_1 = dip(d1_max,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)    
    dip_max = np.max(np.asarray([dip1_0, dip1_1, dip2_0, dip2_1]))
    dip_min = np.min(np.asarray([dip1_0, dip1_1, dip2_0, dip2_1]))
    dip_min = min(dip_min, 0)
    dip_range = 0.05*(dip_max - dip_min)
    
    if drug_id == "drug2":
        left_drug = mcmc_data['drug1_name'].values[0]
        right_drug = drug    
    else:
        left_drug = drug
        right_drug = mcmc_data['drug2_name'].values[0]
    
    ax = fig.add_subplot(111)
    ax.plot(y, dip1_0, c="#ff7f0e", linewidth=lw, label="0uM %s"%right_drug)
    ax.plot(y, dip1_1, c="#2ca02c", linewidth=lw, label="%duM %s"%(np.power(10.,d2_max)*10**6, right_drug))
    ax.plot(y, 0*y, 'k--', alpha=0.7)
    ax.set_ylim(dip_min-dip_range, dip_max+dip_range)
    ax.set_xlim(y[0],y[-1])
    f = interpolate.interp1d(dip1_0,y)
    y1 = (max(dip1_0)-(max(dip1_0)-min(dip1_0))/2)
    x1 = f(y1)    
    plt.scatter(x1,y1,s=20,marker= 'o',color='#ff7f0e',zorder=100)
    plt.plot((x1,x1),(ax.get_ylim()),color = '#ff7f0e',linestyle='--')
    f = interpolate.interp1d(dip1_1,y)
    y1 = (max(dip1_1)-(max(dip1_1)-min(dip1_1))/2)
    x1 = f(y1)    
    plt.scatter(x1,y1,s=20,marker='o',color='#2ca02c',zorder=100)
    plt.plot((x1,x1),(ax.get_ylim()),color = '#2ca02c',linestyle='--')
    ax.set_ylabel(r'DIP rate ($h^{-1}$)')
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.55),
          fancybox=False, shadow=False, ncol=1,handletextpad=0.01)
    ax.set_xlabel(r'log(%s)[M]'%left_drug)

    if fname is None: show()
    else: plt.savefig(right_drug+'_slices_'+fname)
    
    fig = figure(figsize=(1.5,1))
    ax = fig.add_subplot(111)
    ax.plot(t, dip2_0, c="#d62728", linewidth=lw, label="0uM %s"%left_drug)
    ax.plot(t, dip2_1, c="#9467bd", linewidth=lw, label="%duM %s"%(np.power(10.,d1_max)*10**6, left_drug))
    ax.plot(t, 0*t, 'k--', alpha=0.7)
    ax.set_ylim(dip_min-dip_range, dip_max+dip_range)
    ax.set_xlim(t[0],t[-1])
    
    f = interpolate.interp1d(dip2_0,t)
    y1 = (max(dip2_0)-(max(dip2_0)-min(dip2_0))/2)
    x1 = f(y1)    
    plt.scatter(x1,y1,s=20,marker= 'o',color='#d62728',zorder=100)
    plt.plot((x1,x1),(ax.get_ylim()),color = '#d62728',linestyle='--')
    f = interpolate.interp1d(dip2_1,t)
    y1 = (max(dip2_1)-(max(dip2_1)-min(dip2_1))/2)
    x1 = f(y1)    
    plt.scatter(x1,y1,s=20,marker='o',color='#9467bd',zorder=100)
    plt.plot((x1,x1),(ax.get_ylim()),color = '#9467bd',linestyle='--')
    
    ax.set_xlabel("log(%s)[M]"%right_drug)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.55),
          fancybox=False, shadow=False, ncol=1,handletextpad=0.01)
  #    
    if fname is None: show()
    else: plt.savefig(left_drug+'_slices_'+fname)

    
    return dip1_0, dip1_1, dip2_0, dip2_1





############################################################
#Function to plot the interactive surface plots in plotly    
############################################################

# 1D Hill function for fitting
def Edrug1D(d,E0,Em,C,h):
    return Em + (E0-Em) / (1 + (d/C)**h)

# 2D Hill function for fitting
def Edrug2D_DB(d,E0,E1,E2,E3,C1,C2,h1,h2,alpha):
    d1 = d[0]
    d2 = d[1]
    Ed = ( C1**h1*C2**h2*E0 + d1**h1*C2**h2*E1 + C1**h1*d2**h2*E2 + alpha**h2*d1**h1*d2**h2*E3 ) / \
                ( C1**h1*C2**h2 + d1**h1*C2**h2 + C1**h1*d2**h2 + alpha**h2*d1**h1*d2**h2 )
    return Ed
    
#2D dose response surface function not assuming detail balance but assuming proliferation 
# is << than the rate of transition between the different states.
def Edrug2D_NDB(d,E0,E1,E2,E3,r1,r2,C1,C2,h1,h2,alpha1,alpha2):
    d1 = d[0]
    d2 = d[1]
    Ed = (E3*r1*(alpha1*d2)**h2*(alpha2*d1**2)**h1 + 
          E3*r2*(alpha2*d1)**h1*(alpha1*d2**2)**h2 + 
          C2**h2*E0*r1*(C1*alpha2*d1)**h1 + 
          C1**h1*E0*r2*(C2*alpha1*d2)**h2 + 
          C2**h2*E1*r1*(alpha2*d1**2)**h1 + 
          C1**h1*E2*r2*(alpha1*d2**2)**h2 + 
          E3*d2**h2*r1*(C1*alpha2*d1)**h1 + 
          E3*d1**h1*r2*(C2*alpha1*d2)**h2 + 
          C1**(2*h1)*C2**h2*E0*r1 + 
          C1**h1*C2**(2*h2)*E0*r2 + 
          C1**(2*h1)*E2*d2**h2*r1 + 
          C2**(2*h2)*E1*d1**h1*r2 + 
          E1*r2*(C2*d2)**h2*(alpha2*d1)**h1 + 
          E2*r1*(C1*d1)**h1*(alpha1*d2)**h2 +
          C2**h2*E1*r1*(C1*d1)**h1 +
          C1**h1*E2*r2*(C2*d2)**h2)/  \
         (r1*(alpha1*d2)**h2*(alpha2*d1**2)**h1 + 
          r2*(alpha2*d1)**h1*(alpha1*d2**2)**h2 + 
          C2**h2*r1*(C1*alpha2*d1)**h1 +
          C1**h1*r2*(C2*alpha1*d2)**h2 + 
          C2**h2*r1*(alpha2*d1**2)**h1 + 
          C1**h1*r2*(alpha1*d2**2)**h2 + 
          d2**h2*r1*(C1*alpha2*d1)**h1 + 
          d1**h1*r2*(C2*alpha1*d2)**h2 + 
          C1**(2*h1)*C2**h2*r1 + 
          C1**h1*C2**(2*h2)*r2 + 
          C1**(2*h1)*d2**h2*r1 + 
          C2**(2*h2)*d1**h1*r2 + 
          r1*(C1*d1)**h1*(alpha1*d2)**h2 + 
          r2*(C2*d2)**h2*(alpha2*d1)**h1 + 
          C2**h2*r1*(C1*d1)**h1 + 
          C1**h1*r2*(C2*d2)**h2)
         
    return Ed





def DosePlots_PLY(d1,d2,dip,dip_sd,drug1_name,drug2_name,popt3,which_model,target_cell_line,expt,metric_name):
    #To plot doses on log space
    t_d1 = d1.copy()
    t_d2 = d2.copy()
    t_d1[t_d1==0] = min(t_d1[t_d1!=0])/10
    t_d2[t_d2==0] = min(t_d2[t_d2!=0])/10

    trace1 = go.Scatter3d(x=np.log10(t_d1),y=np.log10(t_d2),z=dip,mode='markers',
                        marker=dict(size=3,color = dip,colorscale = 'coolwarm',
                                    line=dict(
                                              color='rgb(0,0,0)',width=1)),
                        )    
    #Plot the fit
    conc1 = 10**np.linspace(np.log10(min(t_d1)),np.log10(max(t_d1)), 30)
    conc2 = 10**np.linspace(np.log10(min(t_d2)), np.log10(max(t_d2)), 30)
    conc3 = np.array([(c1,c2) for c1 in conc1 for c2 in conc2])
    if which_model == 'DB':
        twoD_fit = np.array([Edrug2D_DB(d,*popt3) for d in conc3])
    elif which_model=='NDB':
        twoD_fit = np.array([Edrug2D_NDB(d,*popt3) for d in conc3])
        
    zl = [min(twoD_fit),max(twoD_fit)]
    twoD_fit = np.resize(twoD_fit, (len(conc1),len(conc2)))
    
    [X,Y] = np.meshgrid(np.log10(conc1), np.log10(conc2))

    trace2 = go.Surface(x=X,y=Y,z=twoD_fit.transpose(),
                      colorscale = 'coolwarm',
                      opacity=.8,
                      contours  = dict( 
                                        y = dict(highlight=False,show=False,color='#444',width=1),
                                        x = dict(highlight=False,show=False,color='#444',width = 1 ),
                                        z = dict(highlight=True,show=True,highlightwidth=4,width=3,usecolormap=False)
                                        ))

                      
    
    layout = go.Layout(
                       scene = dict(
                                    xaxis=dict(    
                                               range = [np.log10(min(t_d1)),np.log10(max(t_d1))],
                                               title = 'log(' + drug1_name + ')',
                                               autorange=True,
                                               showgrid=True,
                                               zeroline=False,
                                               showline=False,
                                               ticks='',
                                               showticklabels=True,
                                               gridwidth = 5,
                                               gridcolor = 'k'
                                    ),
                                    yaxis=dict(    
                                               range = [np.log10(min(t_d1)),np.log10(max(t_d1))],
                                               title = 'log(' + drug2_name + ')',
                                               autorange=True,
                                               showgrid=True,
                                               zeroline=False,
                                               showline=False,
                                               ticks='',
                                               showticklabels=True,
                                               gridwidth = 5,
                                               gridcolor = 'k'
                                    ),
                                    zaxis=dict(
                                               range = zl,
                                               title = metric_name,
                                               autorange=False,
                                               tick0=np.min(twoD_fit),
                                               dtick=(np.max(twoD_fit)-np.min(twoD_fit))/5,
                                               showgrid=True,
                                               zeroline=True,
                                               zerolinewidth=5,
                                               showline=False,
                                               ticks='',
                                               showticklabels=True,
                                               gridwidth=5,
                                               gridcolor = 'k'
                                    )
                                    ),
                    margin=dict(
                                l=0,
                                r=0,
                                b=0,
                                t=0
                    ),
                    showlegend=False,
                    font=dict(family='Arial', size=18),
                    title = target_cell_line
    )
        
                    
    data = [trace1,trace2]


    for e,i in enumerate(dip_sd):
        x = np.log10(t_d1[e])*np.ones((2,))
        y = np.log10(t_d2[e])*np.ones((2,))
        z = np.array([dip[e]+i,dip[e]-i])
        trace = go.Scatter3d(
            x=x, y=y, z=z,
            mode = 'lines',
            line=dict(
                color='#1f77b4',
                width=2
            )
        )
        data.append(trace)
    
    fig = go.Figure(data=data,layout=layout)
    
    camera = dict(
                  up=dict(x=0,y=0,z=1),
                  center=dict(x=0,y=0,z=0),
                  eye = dict(x=1.25,y=-1.25,z=1.25))
    fig['layout'].update(scene=dict(camera=camera))
    s_idx = expt.rfind(os.sep)
    expt = expt[s_idx+1:-4]
    plot(fig,filename = 'html/{}_{}_{}_{}_{}_plotly_doseResponseSurface.html'.format(target_cell_line,drug1_name,drug2_name,which_model,expt),auto_open=False)










