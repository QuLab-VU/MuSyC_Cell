# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 16:13:54 2018

@author: xnmeyer
"""
#Import libraries and configure default font/axes parameters
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#Adjust text from  https://github.com/Phlya/adjustText
#Can be installed with pip
import scipy.stats as st
font = {'family' : 'normal',
        'weight':'normal',
        'size'   : 8}
axes = {'linewidth': 2}
rc('font', **font)
rc('axes',**axes)
from scipy import integrate


T = pd.read_csv('../../../Data/MasterResults.csv')
T = T[T['cell_line']!='PC9C1']
T = T[T['cell_line']!='BR1']
T = T[T['cell_line']!='PC9AZR']
T = T[T['MCMC_converge']==1]
T = T[T['model_level']>4]
T = T.reset_index(drop = True)

#marker for each drug combo
#Color for each cell line 
#Markers for each sub class
#Colors for each super class
cell_lines = list(T['cell_line'].unique()) 
num_classes = len(T['cell_line'].unique())
a = []
for letter in range(97,97+16):
    a.append(chr(letter))
markers = []
for i in range(16):
    markers.append('$' + a[i] + '$')    
t = T.drop_duplicates(['drug1_name','drug2_name'])[['drug1_name','drug2_name']]
drugs = list(np.sort(t['drug1_name'] + '_' + t['drug2_name']))
col = ['xkcd:red purple','tab:blue', 'tab:orange', 'tab:red', 'tab:purple','tab:green','xkcd:blue green','tab:brown']
cmaps = ['RdPu','Blues','Oranges','Reds','Purples','Greens','GnBu','YlOrBr']
T['drug_label'] = T['drug1_name']+'_'+T['drug2_name']

################################################################################
###Make initial density plot of the all the drugs
#################################################################################
fig = plt.figure(figsize=(1.2,1.2))
ax = [];
ax.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=3))
ax.append(plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=3))
ax.append(plt.subplot2grid((4,4),(1,3),rowspan=3,colspan=1))
plt.subplots_adjust(wspace=0,hspace=0)
#Alpha values (RAF potentiates MEK)
x = np.array(T['log_alpha1'])
#Beta values
y = np.array(T['beta_obs_norm'])
xx,yy=np.mgrid[min(x):max(x):100j,min(y):max(y):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([x, y])
#Make a 2D gaussain
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
#PLot the 2D guassian contours
tiers = [.01,.05,1]
cfset = ax[0].contourf(xx, yy, f,tiers, cmap='Greys')
cset = ax[0].contour(xx, yy, f,tiers, colors='k')
ax[0].clabel(cset, inline=1, fontsize=8)
#Format axes
xmin,xmax = min(x)-.5,max(x)+.5
ymin,ymax = min(y)-.3,max(y)+.1
ax[0].set_xlim((xmin,xmax))
ax[0].set_ylim((ymin,ymax))
ax[0].set_xlabel(r'$log(\alpha_1)$(RAF potentiates MEK)')
ax[0].set_ylabel(r'$\beta_{obs}$')
ax[0].plot([0,0],(ymin,ymax),linestyle='--',color='k')
ax[0].plot((xmin,xmax),[0,0],linestyle='--',color='k')

#Now plot 1D pdf above and on the side axes corresponding to alpha and beta respectively
#Format axes
ax[1].xaxis.tick_top()
ax[1].get_yaxis().set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_visible(False)
ax[2].yaxis.tick_right()
ax[2].get_xaxis().set_visible(False)
ax[2].spines['top'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[2].spines['left'].set_visible(False)

#Plot alpha pdf
kde = st.kde.gaussian_kde(x)
allD_ax = np.linspace(xmin,xmax,100)
allD_ay = kde(allD_ax)
ax[1].plot(allD_ax,allD_ay,'k',zorder=10)
tiers = [.5,.95]
tier_val_a = []
cnt = .001
for t in tiers:
    while True:
        cnt+=.001
        if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
            tier_val_a.append(cnt)
            break    
ax[1].set_xticks(np.array(tier_val_a))
ax[1].set_xticklabels([str(int(i*100))+'%' for i in tiers])

#Plot beta pdf
kde = st.kde.gaussian_kde(y)
allD_bx = np.linspace(ymin,ymax,100)
allD_by = kde(allD_bx)
ax[2].plot(allD_by,allD_bx,'k',zorder=10)
tiers = [.5,.95]
tier_val_b= []
cnt = .001
for t in tiers:
    while True:
        cnt+=.001
        if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
            tier_val_b.append(cnt)
            break    
ax[2].set_yticks(np.array(tier_val_b))
ax[2].set_yticklabels([str(int(i*100))+'%' for i in tiers])
plt.savefig('Plots/SuppMel_density_All_RAFi.pdf')


################################################################################
###Make density plot for each super and sub classes
#################################################################################
#Density maps for each cell line
for i in range(num_classes):
    fig = plt.figure(figsize=(1.2,1.2))
    ax = [];
    ax.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=3))
    ax.append(plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=3))
    ax.append(plt.subplot2grid((4,4),(1,3),rowspan=3,colspan=1))
    plt.subplots_adjust(wspace=0,hspace=0)
    a1 =  x[np.in1d(T['cell_line'],cell_lines[i])]   
    a2 =  y[np.in1d(T['cell_line'],cell_lines[i])]   
    values = np.vstack([a1, a2])
    #Make a 2D gaussain
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    tiers = [.01,.05,1]
    cfset = ax[0].contourf(xx, yy, f,tiers, cmap=cmaps[i])
    cset = ax[0].contour(xx, yy, f,tiers, colors='k')
    ax[0].clabel(cset, inline=1, fontsize=6)
    #Format axes
    ax[0].set_xlim((xmin,xmax))
    ax[0].set_ylim((ymin,ymax))
    ax[0].plot([0,0],(ymin,ymax),linestyle='--',color='k')
    ax[0].plot((xmin,xmax),[0,0],linestyle='--',color='k')
    ax[1].xaxis.tick_top()
    ax[1].get_yaxis().set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].spines['left'].set_visible(False)    
    ax[2].yaxis.tick_right()
    ax[2].get_xaxis().set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[2].spines['bottom'].set_visible(False)
    ax[2].spines['left'].set_visible(False)
    kde = st.kde.gaussian_kde(a1)
    subD_ax = np.linspace(xmin,xmax,100)
    subD_ay = kde(subD_ax)
    ax[1].plot(subD_ax,subD_ay,col[i],zorder=20)
    tiers = [.5,.95]
    tier_val_a = []
    cnt = .001
    for t in tiers:
        while True:
            cnt+=.001
            if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
                tier_val_a.append(cnt)
                break    
    ax[1].set_xticks(np.array(tier_val_a))
    ax[1].set_xticklabels([str(int(j*100))+'%' for j in tiers])
    ax[1].plot(allD_ax,allD_ay,'k',zorder=10)        
    kde = st.kde.gaussian_kde(a2)
    subD_bx = np.linspace(ymin,ymax,100)
    subD_by = kde(subD_bx)
    ax[2].plot(subD_by,subD_bx,col[i],zorder=20)
    tiers = [.5,.95]
    tier_val_b= []
    cnt = .001
    for t in tiers:
        while True:
            cnt+=.001
            if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
                tier_val_b.append(cnt)
                break    
    ax[2].set_yticks(np.array(tier_val_b))
    ax[2].set_yticklabels([str(int(j*100))+'%' for j in tiers])
    ax[2].plot(allD_by,allD_bx,'k',zorder=10)
    plt.savefig('Plots/SuppDSDs_Melanoma_density_' + cell_lines[i]+'RAFi.pdf') 

############################################################################
#Plot drug synergy diagrams for all, as well as each super and sub classes
############################################################################
#DSDs...
fig = plt.figure(figsize=(8.5,2))
ax = []
#num_classes+1 plots to include all subplot
for i in range(4):
    ax.append(plt.subplot2grid((2,4),(0,i)))
for i in range(4):
    ax.append(plt.subplot2grid((2,4),(1,i)))
plt.subplots_adjust(wspace=0,hspace=.4)
#Empty list of lists to hold the text labels 
text = [ [] for i in range(num_classes)]
for i in range(num_classes):
    plt.sca(ax[i])
    seen = []
    for ind in T[T['cell_line']==cell_lines[i]].index:
        mk = T['drug_label'].loc[ind]
        #If this type of marker has not been seen before, label point for 
        #legend and add plot title
        if ~np.in1d(mk,seen)[0]:
            seen.append(mk)
            ax[i].scatter(x[ind],y[ind],
                        s=30,color=col[i],marker=markers[drugs.index(mk)],
                        zorder=10,edgecolor='k',linewidth=1,label=mk)
            plt.title(cell_lines[i])
        #If the marker has been seen, dont label...
        else:
            ax[i].scatter(x[ind],y[ind],
                        s=30,color=col[i],marker=markers[drugs.index(mk)],
                        zorder=10,edgecolor='k',linewidth=1)
#Format all the axes
for e in range(num_classes):
    ax[e].set_xlim((xmin,xmax))
    ax[e].set_ylim((ymin,ymax))
    y_lim = (ymin,ymax)
    x_lim = (xmin,xmax)
    plt.sca(ax[e])
    plt.plot([0,0],y_lim,linestyle='--',color='k')
    plt.plot(x_lim,[0,0],linestyle='--',color='k')
    plt.sca(ax[e])
    if e==4:
        plt.xlabel(r'$log(\alpha_1)$ (RAF potentiates MEK)')
        plt.ylabel(r'$\beta_{obs}$')
    else:
        ax[e].set_xticklabels([]) 
        ax[e].set_yticklabels([])
    ax[e].tick_params(direction='in')
    
#Adjust text function from 
#pip install adjustText package:  https://github.com/Phlya/adjustText
plt.savefig('Plots/SuppDSDs_Melanoma_CellLine_RAFi.pdf')







##################################################################################
####Make boxplot of E1,E2,E3 for each super class and sub class
##################################################################################
fig = plt.figure(figsize=(8.5,4.5))
ax = []
for i in range(4):
    ax.append(plt.subplot2grid((2,4),(0,i)))
for i in range(4):
    ax.append(plt.subplot2grid((2,4),(1,i)))
plt.subplots_adjust(wspace=0,hspace=.4)

e1 = T['E1_obs'].values
e2 = T['E2_obs'].values                
e3 = T['E3_obs'].values #Combined max effect

#Make plot for all drugs
plt.sca(ax[0])
x1 = 0+(np.random.rand(len(e1))-.5)/5
x2 = 1+(np.random.rand(len(e2))-.5)/5
x3 = 2+(np.random.rand(len(e3))-.5)/5

#make plot for each super class only plotting max drug X (E2) and max drug X+max drug osi (E3)
for i in range(num_classes):
    plt.sca(ax[i])
    seen = []
    indx = T[T['cell_line']==cell_lines[i]].index
    a1 = [e1[j] for j in indx]
    a2 = [e2[j] for j in indx]
    a3 = [e3[j] for j in indx]
    x1 = 0+(np.random.rand(len(a1))-.5)/5
    x2 = 1+(np.random.rand(len(a2))-.5)/5
    x3 = 2+(np.random.rand(len(a3))-.5)/5
    text = []
#    for en,ind in enumerate(indx):
#        plt.scatter(x1[en],e1[ind],
#                    s=30,c=col[i],marker=markers[drugs.index(T['drug_label'].loc[ind])],
#                    zorder=10,edgecolor='k',linewidth=1)
#        plt.scatter(x2[en],e2[ind],
#                    s=30,c=col[i],marker=markers[drugs.index(T['drug_label'].loc[ind])],
#                    zorder=10,edgecolor='k',linewidth=1)
#        plt.scatter(x3[en],e3[ind],
#                    s=30,c=col[i],marker=markers[drugs.index(T['drug_label'].loc[ind])],
#                    zorder=10,edgecolor='k',linewidth=1)
    plt.title(cell_lines[i])
#    for j in range(len(a2)):
#        plt.plot((x2[j],x3[j]),(a2[j],a3[j]),linewidth=.5,color='k')
#        plt.plot((x1[j],x3[j]),(a1[j],a3[j]),linewidth=.5,color='k')
    
    bp = plt.boxplot([a1,a2,a3],positions=[0,1,2],showfliers=False,zorder=10000,boxprops = {'linewidth':2},whiskerprops={'linewidth':2},medianprops = {'linewidth':2})
    ax[i].set_xlim([-.25,2.25])            
    ax[i].set_ylim((min(e3),T['E0'].max()))
    ax[i].set_xticks([0,1,2])
    plt.plot([-.25,2.25],[0,0],linestyle='--',color='k')
    if i==0:
        ax[i].set_xticklabels(['Max RAFi\n(E1)','Max MEKi\n(E2)','RAFi+MEKi\n(E3)'])
        ax[i].set_ylabel(r'DIP $s^{-1}$')
    else:
        ax[i].set_xticklabels([])
        ax[i].set_yticks([])

plt.savefig('Plots/Supp_MelFig_effect.pdf')

###############################################################################
###############################################################################
##############################################################################
#Same thing except with Log alpha2






################################################################################
###Make initial density plot of the all the drugs
#################################################################################
fig = plt.figure(figsize=(1.2,1.2))
ax = [];
ax.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=3))
ax.append(plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=3))
ax.append(plt.subplot2grid((4,4),(1,3),rowspan=3,colspan=1))
plt.subplots_adjust(wspace=0,hspace=0)
#Alpha values (RAF potentiates MEK)
x = np.array(T['log_alpha2'])
#Beta values
y = np.array(T['beta_obs_norm'])
xx,yy=np.mgrid[min(x):max(x):100j,min(y):max(y):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([x, y])
#Make a 2D gaussain
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
#PLot the 2D guassian contours
tiers = [.01,.05,1]
cfset = ax[0].contourf(xx, yy, f,tiers, cmap='Greys')
cset = ax[0].contour(xx, yy, f,tiers, colors='k')
ax[0].clabel(cset, inline=1, fontsize=8)
#Format axes
xmin,xmax = min(x)-.5,max(x)+.5
ymin,ymax = min(y)-.3,max(y)+.1
ax[0].set_xlim((xmin,xmax))
ax[0].set_ylim((ymin,ymax))
ax[0].set_xlabel(r'$log(\alpha_1)$(MEK potentiates RAF)')
ax[0].set_ylabel(r'$\beta_{obs}$')
ax[0].plot([0,0],(ymin,ymax),linestyle='--',color='k')
ax[0].plot((xmin,xmax),[0,0],linestyle='--',color='k')

#Now plot 1D pdf above and on the side axes corresponding to alpha and beta respectively
#Format axes
ax[1].xaxis.tick_top()
ax[1].get_yaxis().set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_visible(False)
ax[2].yaxis.tick_right()
ax[2].get_xaxis().set_visible(False)
ax[2].spines['top'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[2].spines['left'].set_visible(False)

#Plot alpha pdf
kde = st.kde.gaussian_kde(x)
allD_ax = np.linspace(xmin,xmax,100)
allD_ay = kde(allD_ax)
ax[1].plot(allD_ax,allD_ay,'k',zorder=10)
tiers = [.5,.95]
tier_val_a = []
cnt = .001
for t in tiers:
    while True:
        cnt+=.001
        if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
            tier_val_a.append(cnt)
            break    
ax[1].set_xticks(np.array(tier_val_a))
ax[1].set_xticklabels([str(int(i*100))+'%' for i in tiers])

#Plot beta pdf
kde = st.kde.gaussian_kde(y)
allD_bx = np.linspace(ymin,ymax,100)
allD_by = kde(allD_bx)
ax[2].plot(allD_by,allD_bx,'k',zorder=10)
tiers = [.5,.95]
tier_val_b= []
cnt = .001
for t in tiers:
    while True:
        cnt+=.001
        if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
            tier_val_b.append(cnt)
            break    
ax[2].set_yticks(np.array(tier_val_b))
ax[2].set_yticklabels([str(int(i*100))+'%' for i in tiers])
plt.savefig('Plots/SuppMel_density_All_MEKi.pdf')


################################################################################
###Make density plot for each super and sub classes
#################################################################################
#Density maps for each cell line
for i in range(num_classes):
    fig = plt.figure(figsize=(1.2,1.2))
    ax = [];
    ax.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=3))
    ax.append(plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=3))
    ax.append(plt.subplot2grid((4,4),(1,3),rowspan=3,colspan=1))
    plt.subplots_adjust(wspace=0,hspace=0)
    a1 =  x[np.in1d(T['cell_line'],cell_lines[i])]   
    a2 =  y[np.in1d(T['cell_line'],cell_lines[i])]   
    values = np.vstack([a1, a2])
    #Make a 2D gaussain
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    tiers = [.01,.05,1]
    cfset = ax[0].contourf(xx, yy, f,tiers, cmap=cmaps[i])
    cset = ax[0].contour(xx, yy, f,tiers, colors='k')
    ax[0].clabel(cset, inline=1, fontsize=6)
    #Format axes
    ax[0].set_xlim((xmin,xmax))
    ax[0].set_ylim((ymin,ymax))
    ax[0].plot([0,0],(ymin,ymax),linestyle='--',color='k')
    ax[0].plot((xmin,xmax),[0,0],linestyle='--',color='k')
    ax[1].xaxis.tick_top()
    ax[1].get_yaxis().set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].spines['left'].set_visible(False)    
    ax[2].yaxis.tick_right()
    ax[2].get_xaxis().set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[2].spines['bottom'].set_visible(False)
    ax[2].spines['left'].set_visible(False)
    kde = st.kde.gaussian_kde(a1)
    subD_ax = np.linspace(xmin,xmax,100)
    subD_ay = kde(subD_ax)
    ax[1].plot(subD_ax,subD_ay,col[i],zorder=20)
    tiers = [.5,.95]
    tier_val_a = []
    cnt = .001
    for t in tiers:
        while True:
            cnt+=.001
            if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
                tier_val_a.append(cnt)
                break    
    ax[1].set_xticks(np.array(tier_val_a))
    ax[1].set_xticklabels([str(int(j*100))+'%' for j in tiers])
    ax[1].plot(allD_ax,allD_ay,'k',zorder=10)        
    kde = st.kde.gaussian_kde(a2)
    subD_bx = np.linspace(ymin,ymax,100)
    subD_by = kde(subD_bx)
    ax[2].plot(subD_by,subD_bx,col[i],zorder=20)
    tiers = [.5,.95]
    tier_val_b= []
    cnt = .001
    for t in tiers:
        while True:
            cnt+=.001
            if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
                tier_val_b.append(cnt)
                break    
    ax[2].set_yticks(np.array(tier_val_b))
    ax[2].set_yticklabels([str(int(j*100))+'%' for j in tiers])
    ax[2].plot(allD_by,allD_bx,'k',zorder=10)
    plt.savefig('Plots/SuppDSDs_Melanoma_density_' + cell_lines[i]+'MEKi.pdf') 

############################################################################
#Plot drug synergy diagrams for all, as well as each super and sub classes
############################################################################
#DSDs...
fig = plt.figure(figsize=(8.5,2))
ax = []
#num_classes+1 plots to include all subplot
for i in range(4):
    ax.append(plt.subplot2grid((2,4),(0,i)))
for i in range(4):
    ax.append(plt.subplot2grid((2,4),(1,i)))
plt.subplots_adjust(wspace=0,hspace=.4)
#Empty list of lists to hold the text labels 
text = [ [] for i in range(num_classes)]
for i in range(num_classes):
    plt.sca(ax[i])
    seen = []
    for ind in T[T['cell_line']==cell_lines[i]].index:
        mk = T['drug_label'].loc[ind]
        #If this type of marker has not been seen before, label point for 
        #legend and add plot title
        if ~np.in1d(mk,seen)[0]:
            seen.append(mk)
            ax[i].scatter(x[ind],y[ind],
                        s=30,color=col[i],marker=markers[drugs.index(mk)],
                        zorder=10,edgecolor='k',linewidth=1,label=mk)
            plt.title(cell_lines[i])
        #If the marker has been seen, dont label...
        else:
            ax[i].scatter(x[ind],y[ind],
                        s=30,color=col[i],marker=markers[drugs.index(mk)],
                        zorder=10,edgecolor='k',linewidth=1)
#Format all the axes
for e in range(num_classes):
    ax[e].set_xlim((xmin,xmax))
    ax[e].set_ylim((ymin,ymax))
    y_lim = (ymin,ymax)
    x_lim = (xmin,xmax)
    plt.sca(ax[e])
    plt.plot([0,0],y_lim,linestyle='--',color='k')
    plt.plot(x_lim,[0,0],linestyle='--',color='k')
    plt.sca(ax[e])
    if e==4:
        plt.xlabel(r'$log(\alpha_1)$ (MEK potentiates RAF)')
        plt.ylabel(r'$\beta_{obs}$')
    else:
        ax[e].set_xticklabels([]) 
        ax[e].set_yticklabels([])
    ax[e].tick_params(direction='in')
    
#Adjust text function from 
#pip install adjustText package:  https://github.com/Phlya/adjustText
plt.savefig('Plots/SuppDSDs_Melanoma_CellLine_MEKi.pdf')










#####################################################################
#Do alpha alpha plots for cell-lines
######

################################################################################
###Make initial density plot of the all the drugs
#################################################################################
fig = plt.figure(figsize=(1.2,1.2))
ax = [];
ax.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=3))
ax.append(plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=3))
ax.append(plt.subplot2grid((4,4),(1,3),rowspan=3,colspan=1))
plt.subplots_adjust(wspace=0,hspace=0)
#Alpha values (RAF potentiates MEK)
x = np.array(T['log_alpha1'])
#Beta values
y = np.array(T['log_alpha2'])
xx,yy=np.mgrid[min(x):max(x):100j,min(y):max(y):100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([x, y])
#Make a 2D gaussain
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)
#PLot the 2D guassian contours
tiers = [.001,.01,1]
cfset = ax[0].contourf(xx, yy, f,tiers, cmap='Greys')
cset = ax[0].contour(xx, yy, f,tiers, colors='k')
ax[0].clabel(cset, inline=1, fontsize=8)
#Format axes
xmin,xmax = min(x)-.5,max(x)+.5
ymin,ymax = min(y)-.3,max(y)+.1
ax[0].set_xlim((xmin,xmax))
ax[0].set_ylim((ymin,ymax))
ax[0].set_xlabel(r'$log(\alpha_1)$(RAF potentiates MEK)')
ax[0].set_ylabel(r'$log(\alpha_2)$(MEK potentiates RAF)')
ax[0].plot([0,0],(ymin,ymax),linestyle='--',color='k')
ax[0].plot((xmin,xmax),[0,0],linestyle='--',color='k')

#Now plot 1D pdf above and on the side axes corresponding to alpha and beta respectively
#Format axes
ax[1].xaxis.tick_top()
ax[1].get_yaxis().set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)
ax[1].spines['left'].set_visible(False)
ax[2].yaxis.tick_right()
ax[2].get_xaxis().set_visible(False)
ax[2].spines['top'].set_visible(False)
ax[2].spines['bottom'].set_visible(False)
ax[2].spines['left'].set_visible(False)

#Plot alpha pdf
kde = st.kde.gaussian_kde(x)
allD_ax = np.linspace(xmin,xmax,100)
allD_ay = kde(allD_ax)
ax[1].plot(allD_ax,allD_ay,'k',zorder=10)
tiers = [.5,.95]
tier_val_a = []
cnt = .001
for t in tiers:
    while True:
        cnt+=.001
        if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
            tier_val_a.append(cnt)
            break    
ax[1].set_xticks(np.array(tier_val_a))
ax[1].set_xticklabels([str(int(i*100))+'%' for i in tiers])

#Plot beta pdf
kde = st.kde.gaussian_kde(y)
allD_bx = np.linspace(ymin,ymax,100)
allD_by = kde(allD_bx)
ax[2].plot(allD_by,allD_bx,'k',zorder=10)
tiers = [.5,.95]
tier_val_b= []
cnt = .001
for t in tiers:
    while True:
        cnt+=.001
        if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
            tier_val_b.append(cnt)
            break    
ax[2].set_yticks(np.array(tier_val_b))
ax[2].set_yticklabels([str(int(i*100))+'%' for i in tiers])
plt.savefig('Plots/SuppMel_a_a_density_All_RAFivMEKi.pdf')


################################################################################
###Make density plot for each super and sub classes
#################################################################################
#Density maps for each cell line
for i in range(num_classes):
    fig = plt.figure(figsize=(1.2,1.2))
    ax = [];
    ax.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=3))
    ax.append(plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=3))
    ax.append(plt.subplot2grid((4,4),(1,3),rowspan=3,colspan=1))
    plt.subplots_adjust(wspace=0,hspace=0)
    a1 =  x[np.in1d(T['cell_line'],cell_lines[i])]   
    a2 =  y[np.in1d(T['cell_line'],cell_lines[i])]   
    values = np.vstack([a1, a2])
    #Make a 2D gaussain
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    tiers = [.001,.01,1]
    cfset = ax[0].contourf(xx, yy, f,tiers, cmap=cmaps[i])
    cset = ax[0].contour(xx, yy, f,tiers, colors='k')
    ax[0].clabel(cset, inline=1, fontsize=6)
    #Format axes
    ax[0].set_xlim((xmin,xmax))
    ax[0].set_ylim((ymin,ymax))
    ax[0].plot([0,0],(ymin,ymax),linestyle='--',color='k')
    ax[0].plot((xmin,xmax),[0,0],linestyle='--',color='k')
    ax[1].xaxis.tick_top()
    ax[1].get_yaxis().set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].spines['left'].set_visible(False)    
    ax[2].yaxis.tick_right()
    ax[2].get_xaxis().set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[2].spines['bottom'].set_visible(False)
    ax[2].spines['left'].set_visible(False)
    kde = st.kde.gaussian_kde(a1)
    subD_ax = np.linspace(xmin,xmax,100)
    subD_ay = kde(subD_ax)
    ax[1].plot(subD_ax,subD_ay,col[i],zorder=20)
    tiers = [.5,.95]
    tier_val_a = []
    cnt = .001
    for t in tiers:
        while True:
            cnt+=.001
            if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
                tier_val_a.append(cnt)
                break    
    ax[1].set_xticks(np.array(tier_val_a))
    ax[1].set_xticklabels([str(int(j*100))+'%' for j in tiers])
    ax[1].plot(allD_ax,allD_ay,'k',zorder=10)        
    kde = st.kde.gaussian_kde(a2)
    subD_bx = np.linspace(ymin,ymax,100)
    subD_by = kde(subD_bx)
    ax[2].plot(subD_by,subD_bx,col[i],zorder=20)
    tiers = [.5,.95]
    tier_val_b= []
    cnt = .001
    for t in tiers:
        while True:
            cnt+=.001
            if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
                tier_val_b.append(cnt)
                break    
    ax[2].set_yticks(np.array(tier_val_b))
    ax[2].set_yticklabels([str(int(j*100))+'%' for j in tiers])
    ax[2].plot(allD_by,allD_bx,'k',zorder=10)
    plt.savefig('Plots/SuppMel_a_a_densityByCellLine_' + cell_lines[i]+'RAFivMEKi.pdf') 



################################################################################
###Make density plot for each super and sub classes
#################################################################################
#Density maps for each cell line
cmaps = ['RdPu','Blues','Oranges','Reds','Purples','Greens','GnBu','YlOrBr','OrRd','PuRd','BuPu','PuBu','YlGnBu','PuBuGn','BuGn','YlGn']
for i in range(16):
    fig = plt.figure(figsize=(1.2,1.2))
    ax = [];
    ax.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=3))
    ax.append(plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=3))
    ax.append(plt.subplot2grid((4,4),(1,3),rowspan=3,colspan=1))
    plt.subplots_adjust(wspace=0,hspace=0)
    a1 =  x[np.in1d(T['drug_label'],drugs[i])]   
    a2 =  y[np.in1d(T['drug_label'],drugs[i])]   
    values = np.vstack([a1, a2])
    #Make a 2D gaussain
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    tiers = [.001,.01,1]
    cfset = ax[0].contourf(xx, yy, f,tiers, cmap=cmaps[i])
    cset = ax[0].contour(xx, yy, f,tiers, colors='k')
    ax[0].clabel(cset, inline=1, fontsize=6)
    #Format axes
    ax[0].set_xlim((xmin,xmax))
    ax[0].set_ylim((ymin,ymax))
    ax[0].plot([0,0],(ymin,ymax),linestyle='--',color='k')
    ax[0].plot((xmin,xmax),[0,0],linestyle='--',color='k')
    ax[1].xaxis.tick_top()
    ax[1].get_yaxis().set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].spines['left'].set_visible(False)    
    ax[2].yaxis.tick_right()
    ax[2].get_xaxis().set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[2].spines['bottom'].set_visible(False)
    ax[2].spines['left'].set_visible(False)
    kde = st.kde.gaussian_kde(a1)
    subD_ax = np.linspace(xmin,xmax,100)
    subD_ay = kde(subD_ax)
    ax[1].plot(subD_ax,subD_ay,col[i],zorder=20)
    tiers = [.5,.95]
    tier_val_a = []
    cnt = .001
    for t in tiers:
        while True:
            cnt+=.001
            if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
                tier_val_a.append(cnt)
                break    
    ax[1].set_xticks(np.array(tier_val_a))
    ax[1].set_xticklabels([str(int(j*100))+'%' for j in tiers])
    ax[1].plot(allD_ax,allD_ay,'k',zorder=10)        
    kde = st.kde.gaussian_kde(a2)
    subD_bx = np.linspace(ymin,ymax,100)
    subD_by = kde(subD_bx)
    ax[2].plot(subD_by,subD_bx,col[i],zorder=20)
    tiers = [.5,.95]
    tier_val_b= []
    cnt = .001
    for t in tiers:
        while True:
            cnt+=.001
            if integrate.quad(kde,0-cnt,0+cnt)[0]>t:
                tier_val_b.append(cnt)
                break    
    ax[2].set_yticks(np.array(tier_val_b))
    ax[2].set_yticklabels([str(int(j*100))+'%' for j in tiers])
    ax[2].plot(allD_by,allD_bx,'k',zorder=10)
    plt.savefig('Plots/SuppMel_a_a_densityByDrugType_' + drugs[i]+'RAFivMEKi.pdf') 
