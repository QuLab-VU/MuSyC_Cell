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
axes = {'linewidth': 2}
rc('font', **font)
rc('axes',**axes)


##########################################
### D E F I N E   F U N C T I O N S
##########################################    
#incl_gamma=True; data_pts=False; N=50; elev=20; azim=19; fname=None; zlim=(-.1,1.1); zero_conc=1; alpha = 0.9;metric_name="per eff";path=path = 'conflating_plot/Malarial_combinations/';title=None

def plot_surface(fit_params, incl_gamma=False, data_pts=None, N=50, elev=20, azim=19, fname=None, zlim=None, zero_conc=1, alpha = 0.9,metric_name="DIP rate" + r'$h^{-1}$',path=None,title=None):

    e0, e1, e2, e3, h1, h2, r1, r1r, r2, r2r, ec50_1, ec50_2, alpha1, alpha2, gamma1,gamma2,d1_min, d1_max, d2_min, d2_max, drug1_name, drug2_name, expt = get_params(fit_params)
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
        
    # Plot it once on a throwaway axis, to get cbar without alpha problems
    my_cmap_rgb = plt.get_cmap('bwr')(np.arange(256))

    for i in range(3): # Do not include the last column!
        my_cmap_rgb[:,i] = (1 - alpha) + alpha*my_cmap_rgb[:,i]
    my_cmap = ListedColormap(my_cmap_rgb, name='my_cmap')
    surf = ax.plot_surface(yy, tt, zz, cstride=1, rstride=1, cmap = my_cmap, vmin=zmin, vmax=zmax, linewidth=0)
#    cbar_ax = fig.add_axes([0.1, 0.9, 0.6, 0.05])
    cbar = fig.colorbar(surf,ticks=[zmin, 0, zmax],pad=-.08,fraction=.025)
    cbar.ax.set_yticklabels(['%.2f'%zmin, 0, '%.2f'%zmax])
    cbar.solids.set_rasterized(True)
    cbar.solids.set_edgecolor('face')
    cbar.outline.set_visible(False)
    
    plt.cla()
    ax.set_zlim(zmin,zmax)
    # Plot the surface for real, with appropriate alpha
    surf = ax.plot_surface(tt, yy, zz, cstride=1, rstride=1, alpha=alpha, cmap = cm.bwr, vmin=zmin, vmax=zmax, linewidth=0)
    
    # colored curves on left and right
    lw = 5
    if incl_gamma:
        ax.plot(d2_min*np.ones(y.shape)-zero_conc, y, Edrug2D_NDB_hill((np.power(10.,y),np.power(10.,d2_min*np.ones(y.shape)-zero_conc)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), linewidth=lw,color='b')
        ax.plot(d2_max*np.ones(y.shape), y, Edrug2D_NDB_hill((np.power(10.,y),np.power(10.,d2_max*np.ones(y.shape))),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), linewidth=lw,color='b',linestyle=':')

        ax.plot(t, d1_min*np.ones(y.shape)-zero_conc, Edrug2D_NDB_hill((np.power(10.,d1_min*np.ones(t.shape)-zero_conc),np.power(10.,t)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), linewidth=lw,color='r')
        ax.plot(t, d1_max*np.ones(y.shape), Edrug2D_NDB_hill((np.power(10.,d1_max*np.ones(t.shape)),np.power(10.,t)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2), linewidth=lw,color='r',linestyle=':')
    
    else:
        ax.plot(d2_min*np.ones(y.shape)-zero_conc, y, dip(y,d2_min-zero_conc, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), linewidth=lw,color='b')
        ax.plot(d2_max*np.ones(y.shape), y, dip(y,d2_max, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), linewidth=lw,color='b',linestyle=':')
    
        ax.plot(t, d1_min*np.ones(y.shape)-zero_conc, dip(d1_min-zero_conc,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), linewidth=lw,color='r')
        ax.plot(t, d1_max*np.ones(y.shape), dip(d1_max,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r), linewidth=lw,color='r',linestyle=':')
        
    
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
    ax.set_ylabel("log(%s)[M]"%drug1_name,labelpad=-5)   
    ax.set_xlabel("log(%s)[M]"%drug2_name,labelpad=-5)
    ax.set_zlabel(metric_name,labelpad=-3)

    if data_pts is not None:
        # Scatter points, and get error_bars
        scat_d1,scat_d2,scat_dip,scat_erb = read_in_data(fit_params,exp_path=path)
        scat_d1 = np.log10(scat_d1)
        scat_d1[scat_d1==-np.inf] = scat_d1[scat_d1!=-np.inf].min()-zero_conc
        scat_d2 = np.log10(scat_d2)
        scat_d2[scat_d2==-np.inf] = scat_d2[scat_d2!=-np.inf].min()-zero_conc
        ax.scatter(scat_d2, scat_d1, scat_dip, s=1,c='k', depthshade=True)
        # Plot error bars
        for _d1, _d2, _dip, _erb in zip(scat_d1, scat_d2, scat_dip, scat_erb):
            ax.plot([_d2,_d2], [_d1,_d1], [_dip-_erb, _dip+_erb], 'k-', alpha=0.3,linewidth=1)

    #format axis
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    np.linspace(np.ceil(2.3),np.floor(8.4),-np.ceil(2.3)+np.floor(8.4)+1)
    x_ax = np.linspace(np.ceil(d2_min-zero_conc),np.floor(d2_max),int(np.floor(d2_max)-np.ceil(d2_min-zero_conc)+1.))
    x_ax=x_ax[0:-1:np.round(len(x_ax)/3)]
    y_ax = np.linspace(np.ceil(d1_min-zero_conc),np.floor(d1_max),int(np.floor(d1_max)-np.ceil(d1_min-zero_conc)+1.))
    y_ax=y_ax[0:-1:np.round(len(y_ax)/3)]
    ax.set_xticks(x_ax)
    ax.set_yticks(y_ax)
    ax.set_zticks([zmin,0,zmax])
    ax.xaxis.set_tick_params(pad=-4)
    ax.yaxis.set_tick_params(pad=-4)
    ax.zaxis.set_tick_params(pad=-2)


    # Plane of DIP=0
    c_plane = colorConverter.to_rgba('k', alpha=0.3)
    verts = [np.array([(d2_min-zero_conc,d1_min-zero_conc), (d2_min-zero_conc,d1_max), (d2_max,d1_max), (d2_max,d1_min-zero_conc), (d2_min-zero_conc,d1_min-zero_conc)])]
    poly = PolyCollection(verts, facecolors=c_plane)
    ax.add_collection3d(poly, zs=[0], zdir='z')

    # Plot intersection of surface with DIP=0 plane
    plt.contour(tt,yy,zz,levels=[0], linewidths=3, colors='k')

    # If needed, manually set zlim
    if zlim is not None: ax.set_zlim(zlim)
    plt.tight_layout()

    if title is not None:      
        plt.title(title)


    # Save plot, or show it
    if fname is None: plt.show()
    else: plt.savefig(fname+'.pdf',pad_inches=0.,format='pdf');plt.close()



#    #######################
#    #Now plot slices
#    ######################
#    if incl_gamma:
#        dip2_0 = Edrug2D_NDB_hill((np.power(10.,t),np.power(10.,d1_min*np.ones(len(y))-zero_conc)), e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)
#        dip2_1 = Edrug2D_NDB_hill((np.power(10.,t),np.power(10.,d1_max*np.ones(len(y))))          , e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)
#        dip1_0 = Edrug2D_NDB_hill((np.power(10.,d2_min*np.ones(len(y))-zero_conc),np.power(10.,y)), e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)
#        dip1_1 = Edrug2D_NDB_hill((np.power(10.,d2_max*np.ones(len(y))          ),np.power(10.,y)), e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2) 
#    else:
#        dip1_0 = dip(y,d2_min-zero_conc*np.ones(len(y)), e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
#        dip1_1 = dip(y,d2_max*np.ones(len(y)), e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
#        dip2_0 = dip(d1_min-zero_conc*np.ones(len(y)),t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
#        dip2_1 = dip(d1_max*np.ones(len(y)),t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)    
#
#    # colored curves on left and right
#    fig = plt.figure(figsize=(2.5,2.75))
#    lw = 3
#    ax1 = plt.subplot2grid((2,1),(0,0))
#    ax2 = plt.subplot2grid((2,1),(1,0))
#    plt.subplots_adjust(hspace=.1)
#    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
#    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
#    
#    ax1.plot(y, dip1_0, linewidth=lw,color='b',label='0nM '+drug2_name)
#    ax1.plot(y, dip1_1, linewidth=lw,color='b',linestyle=':',label='50nM '+drug2_name)
#    ax1.plot(y,0*y,'k--',alpha=.7)
#    ax1.set_ylim(zmin-.01,zmax+.01)
#    ax1.set_xlim(y[0],y[-1])
#    ax1.set_ylabel(metric_name,labelpad=-2)
#    ax1.set_xlabel("log(%s)[M]"%drug1_name[0:3],labelpad=-1) 
#    ax1.spines['top'].set_visible(False)
#    ax1.spines['right'].set_visible(False)
#    # Put a legend below current axis
#    ax1.legend(loc='lower left',fancybox=True, shadow=False, ncol=1,handletextpad=.5)
#    ax1.set_xticks(x_ax)
#    ax1.set_yticks([zmin,0,zmax])
#    ax1.xaxis.set_tick_params(pad=1)
#    ax1.yaxis.set_tick_params(pad=-1)
#    
#    ax2.plot(t, dip2_0, linewidth=lw,color='r',label='min('+drug1_name[0:3]+')')
#    ax2.plot(t, dip2_1, linewidth=lw,color='r',linestyle=':',label='max('+drug1_name[0:3]+')')
#    ax2.plot(y,0*y,'k--',alpha=.7)
#    ax2.set_ylim(zmin-.01,zmax+.01)
#    ax2.set_xlim(y[0],y[-1])
#    ax2.set_ylabel('',labelpad=-3)
#    ax2.set_xlabel("log(%s)[M]"%drug2_name[0:3],labelpad=-1)
#    ax2.spines['top'].set_visible(False)
#    ax2.spines['right'].set_visible(False)
#    ax2.legend(loc='lower left',fancybox=True, shadow=False, ncol=1,handletextpad=.5)
#    ax2.set_xticks(y_ax)
#    ax2.set_yticks([zmin,0,zmax])
#    ax2.xaxis.set_tick_params(pad=1)
#    ax2.yaxis.set_tick_params(pad=-1)
#    
#    #Plot the EC50
#    plt.sca(ax1)
#    f = interpolate.interp1d(dip1_0,y)
#    y1 = (max(dip1_0)-(max(dip1_0)-min(dip1_0))/2)
#    x1 = f(y1)    
#    plt.scatter(x1,y1,s=30,marker= 'o',edgecolor='b',facecolor='w',zorder=100)
#    plt.plot((x1,x1),(ax1.get_ylim()),color = 'b',linestyle='--')
#    f = interpolate.interp1d(dip1_1,y)
#    y1 = (max(dip1_1)-(max(dip1_1)-min(dip1_1))/2)
#    x1 = f(y1)    
#    plt.scatter(x1,y1,s=20,marker='o',color='b',zorder=100)
#    plt.plot((x1,x1),(ax1.get_ylim()),color = 'b',linestyle='--')
#
#
#    plt.sca(ax2)    
#    f = interpolate.interp1d(dip2_0,t)
#    y1 = (max(dip2_0)-(max(dip2_0)-min(dip2_0))/2)
#    x1 = f(y1)    
#    plt.scatter(x1,y1,s=30,marker= 'o',facecolor='w',edgecolor='r',zorder=100)
#    plt.plot((x1,x1),(ax2.get_ylim()),color = 'r',linestyle='--')
#    f = interpolate.interp1d(dip2_1,t)
#    y1 = (max(dip2_1)-(max(dip2_1)-min(dip2_1))/2)
#    x1 = f(y1)    
#    plt.scatter(x1,y1,s=20,marker='o',color='r',zorder=100)
#    plt.plot((x1,x1),(ax2.get_ylim()),color = 'r',linestyle='--')
#    plt.tight_layout()
#    if fname is None: plt.show()
#    else: plt.savefig(fname+'_slices.pdf',pad_inches=0.,format='pdf');plt.close()

    
    
def plot_slices(fit_params, incl_gamma=False, dd=.1,gN=50, N=10, fname=None, zlim=None, zero_conc=1,metric_name="DIP rate" + r'$h{-1}$',title=None):
    
    e0, e1, e2, e3, h1, h2, r1, r1r, r2, r2r, ec50_1, ec50_2, alpha1, alpha2, gamma1,gamma2,d1_min, d1_max, d2_min, d2_max, drug1_name, drug2_name, expt = get_params(fit_params)
    #ax.set_axis_off()
    y = np.linspace(d1_min-zero_conc, d1_max ,gN)
    t = np.linspace(d2_min-zero_conc, d2_max ,gN)
    
    y_p = np.linspace(d1_min-zero_conc, d1_max ,N)
    t_p = np.linspace(d2_min-zero_conc, d2_max ,N)

    if incl_gamma:
        tt, yy = np.meshgrid(t, y)
        conc3 = np.array([(np.power(10.,c1),np.power(10.,c2)) for c1 in y for c2 in t])
        zz = np.array([Edrug2D_NDB_hill(d,e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2) for d in conc3])
        zz = zz.reshape(gN,gN)
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
        
    # colored curves on left and right
    fig = plt.figure(figsize=(2,2))
    ax1 = plt.subplot2grid((2,10),(0,0),colspan=8)
    ax2 = plt.subplot2grid((2,10),(1,0),colspan=8)
    ax11 = plt.subplot2grid((2,10),(0,9))
    ax22 = plt.subplot2grid((2,10),(1,9))
    ax111 = plt.subplot2grid((2,10),(0,8))
    ax222 = plt.subplot2grid((2,10),(1,8))

    plt.subplots_adjust(hspace=.5,wspace=.7)
    
    cmap1 = cm.Blues(np.linspace(.25,1,N))
    cmap2 = cm.Reds(np.linspace(.25,1,N))
    dx = 1e-6
    if incl_gamma:
        # light grey grid across surface
        for j,ttt in enumerate(t_p):
            dip1_0 = Edrug2D_NDB_hill((np.power(10.,y),np.power(10.,ttt)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)
            f = interpolate.interp1d(dip1_0,y)
            y1 = (max(dip1_0)-(max(dip1_0)-min(dip1_0))/2)
            x1 = f(y1)    
            ax1.scatter(x1,y1,s=20,marker= 'o',color=cmap1[j,:],zorder=100)
            if j == 0 or j==N-1:
                ax1.plot(y, dip1_0, color=cmap1[j,:], linewidth=1)
            #for derivative
            dy =  Edrug2D_NDB_hill((np.power(10.,x1+dx),np.power(10.,ttt)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)-Edrug2D_NDB_hill((np.power(10.,x1-dx),np.power(10.,ttt)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)
            slope = dy/dx
            ax1.plot(np.array([x1-dd,x1+dd]),np.array([y1-dd*slope,y1+dd*slope]),color=cmap1[j,:])
            
        for j,yyy in enumerate(y_p):
            dip1_0 = Edrug2D_NDB_hill((np.power(10.,yyy),np.power(10.,t)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)
            f = interpolate.interp1d(dip1_0,t)
            y1 = (max(dip1_0)-(max(dip1_0)-min(dip1_0))/2)
            x1 = f(y1)    
            ax2.scatter(x1,y1,s=20,marker= 'o',color=cmap2[j,:],zorder=100)
            if j == 0 or j==N-1:
                ax2.plot(t, dip1_0, color=cmap2[j,:], linewidth=1)
            #for derivative
            dy =  Edrug2D_NDB_hill((np.power(10.,yyy),np.power(10.,x1+dx)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)-Edrug2D_NDB_hill((np.power(10.,yyy),np.power(10.,x1-dx)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)
            slope = dy/dx
            ax2.plot(np.array([x1-dd,x1+dd]),np.array([y1-dd*slope,y1+dd*slope]),color=cmap2[j,:])

        
    else:
        for j,ttt in enumerate(t_p):
            dip1_0 = dip(y,ttt, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
            f = interpolate.interp1d(dip1_0,y)
            y1 = (max(dip1_0)-(max(dip1_0)-min(dip1_0))/2)
            x1 = f(y1)    
            ax1.scatter(x1,y1,s=20,marker= 'o',color=cmap1[j,:],zorder=100)
            if j == 0 or j==N-1:
                ax1.plot(y, dip1_0, color=cmap1[j,:], linewidth=1)
            #for derivative
            dy =  Edrug2D_NDB_hill((np.power(10.,x1+dx),np.power(10.,ttt)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)-Edrug2D_NDB_hill((np.power(10.,x1-dx),np.power(10.,ttt)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)
            slope = dy/dx
            ax1.plot(np.array([x1-dd,x1+dd]),np.array([y1-dd*slope,y1+dd*slope]),color=cmap1[j,:])
            
        for j,yyy in enumerate(y_p):
            dip1_0 = dip(yyy,t, e0, e1, e2, e3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
            f = interpolate.interp1d(dip1_0,t)
            y1 = (max(dip1_0)-(max(dip1_0)-min(dip1_0))/2)
            x1 = f(y1)    
            ax2.scatter(x1,y1,s=20,marker= 'o',color=cmap2[j,:],zorder=100)
            if j == 0 or j==N-1:
                ax2.plot(t, dip1_0, color=cmap2[j,:], linewidth=1)
            #for derivative
            dy =  Edrug2D_NDB_hill((np.power(10.,yyy),np.power(10.,t+dx)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)-Edrug2D_NDB_hill((np.power(10.,yyy),np.power(10.,t-dx)),e0,e1,e2,e3,r1,r2,ec50_1,ec50_2,h1,h2,alpha1,alpha2,gamma1,gamma2)
            slope = dy/dx
            ax2.plot(np.array([x1-dd,x1+dd]),np.array([y1-dd*slope,y1+dd*slope]),color=cmap2[j,:])   
    
    
    ax1.set_xlim(d1_min,d1_max)
    ax2.set_xlim(d2_min,d2_max)
    ax1.set_ylim(zmin,zmax)
    ax2.set_ylim(zmin,zmax)    
    ax2.spines['right'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.set_xlabel('log('+drug1_name+')',labelpad=-1)
    ax2.set_xlabel('log('+drug2_name+')',labelpad=-1)
    ax1.set_ylabel(metric_name,labelpad=-1)
    
    ax111.set_xlim(0,1)
    ax111.set_ylim(zmin,zmax)
    ax111.set_xticks([.5])
    ax111.scatter(.5,e2,s=20,c='r',marker='D')
    ax111.spines['right'].set_visible(False)
    ax111.spines['top'].set_visible(False)
    ax111.spines['left'].set_visible(False)
    ax111.set_xticklabels(['max('+drug2_name+')'])
    ax111.set_yticks([])

    ax222.set_xlim(0,1)
    ax222.set_ylim(zmin,zmax)
    ax222.set_xticks([.5])
    ax222.scatter(.5,e1,s=20,c='b',marker='D')
    ax222.spines['right'].set_visible(False)
    ax222.spines['top'].set_visible(False)
    ax222.spines['left'].set_visible(False)
    ax222.set_xticklabels(['max('+drug1_name+')'])
    ax222.set_yticks([])

    cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=mpl.cm.Blues,
                                    norm=mpl.colors.Normalize(vmin=d2_min, vmax=d2_max),
                                    orientation='vertical',drawedges=False)
    
    cb1.set_label('log('+drug2_name+')',labelpad=-1)
    cb1.set_ticks((d2_min,d2_max))
    cb1.outline.set_edgecolor(None)

    cb2 = mpl.colorbar.ColorbarBase(ax22, cmap=mpl.cm.Reds,
                                    norm=mpl.colors.Normalize(vmin=d1_min, vmax=d1_max),
                                    orientation='vertical')
    cb2.set_label('log('+drug1_name+')',labelpad=-1)
    cb2.set_ticks((d1_min,d1_max))    
    cb2.outline.set_edgecolor(None)
    
    ax22.spines['left'].set_visible(False)
    ax11.spines['left'].set_visible(False)
    ax22.spines['top'].set_visible(False)
    ax11.spines['top'].set_visible(False)
    ax22.spines['bottom'].set_visible(False)
    ax11.spines['bottom'].set_visible(False)    

    if fname is None: plt.show()
    else: plt.savefig(fname+'_slices.pdf',pad_inches=0.,format='pdf');plt.close()

    
    
###########################################################################
###Function to read in the data and subset
###########################################################################
def read_in_data(fit_params,exp_path=None):
    expt = exp_path+fit_params['expt']#Experiment file
    target_cell_line = fit_params['cell_line'].str.upper().values[0] #The cell line to consider
    
    #Read in the data into a pandas data frame
    data = pd.read_table(expt.values[0], delimiter=',')
    #format to match the expectation from the fitted data
    data['sample'] = data['sample'].str.upper()
    data['drug1'] = data['drug1'].str.lower()
    data['drug2'] = data['drug2'].str.lower()
    
    #Subset by target cell line
    sub_data = data[data['sample']==target_cell_line]
    #data = data.reset_index()
    
    drug1_name = fit_params['drug1_name'].str.lower().values[0]#Drug1 name to replace in array gen function
    drug2_name = fit_params['drug2_name'].str.lower().values[0]#Drug2 name
    
    #control data for drug 1.  Either
    indx    = ((sub_data['drug1']==drug1_name) & (sub_data['drug2.conc']==0))  & (sub_data['drug1.conc']!=0)
    d1      = sub_data[indx]['drug1.conc'].values
    d2      = sub_data[indx]['drug2.conc'].values
    dip     = sub_data[indx]['effect'].values
    dip_95ci= sub_data[indx]['effect.95ci'].values
    
    indx    = ((sub_data['drug2']==drug1_name) & (sub_data['drug1.conc']==0))  & (sub_data['drug2.conc']!=0)
    d1      = np.concatenate([d1,sub_data[indx]['drug2.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug1.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    #control for drug 2
    indx    = ((sub_data['drug1']==drug2_name) & (sub_data['drug2.conc']==0))  & (sub_data['drug1.conc']!=0)
    d1      = np.concatenate([d1,sub_data[indx]['drug2.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug1.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    indx    = ((sub_data['drug2']==drug2_name) & (sub_data['drug1.conc']==0))  & (sub_data['drug2.conc']!=0)
    d1      = np.concatenate([d1,sub_data[indx]['drug1.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug2.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    #Combination experiment
    indx    = ((sub_data['drug1']==drug1_name) & (sub_data['drug2']==drug2_name)) & ((sub_data['drug1.conc']!=0) & (sub_data['drug2.conc']!=0)) 
    d1      = np.concatenate([d1,sub_data[indx]['drug1.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug2.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    indx    = ((sub_data['drug2']==drug1_name) & (sub_data['drug1']==drug2_name)) & ((sub_data['drug1.conc']!=0) & (sub_data['drug2.conc']!=0)) 
    d1      = np.concatenate([d1,sub_data[indx]['drug2.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug1.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    #The double control condition
    indx    = ((sub_data['drug1.conc']==0) & (sub_data['drug2.conc']==0))
    d1      = np.concatenate([d1,sub_data[indx]['drug1.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug2.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    #Set as standard deviation
    dip_sd = dip_95ci/(2*1.96)
    
    #Remove nan values
    d1      = d1[~np.isnan(dip)]
    d2      = d2[~np.isnan(dip)]
    dip_sd  = dip_sd[~np.isnan(dip)]
    dip     = dip[~np.isnan(dip)]
    
    return d1,d2,dip,dip_sd

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




