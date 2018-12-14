#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Christian T Meyer
email christian.t.meyer@vanderbilt.edu
Code for creating figure 2 in paper 'Quantifying Drug Synergy with respect to potency and efficacy'
Commented and up to date 4-27-18
"""
#Import libraries and configure default font/axes parameters
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

#####################################################################
####Read in and filter the PC9 data...
#####################################################################
T = pd.read_csv('../../Data/MasterResults.csv')
t_pc9 = T[np.logical_or(T['cell_line']=='PC9C1',T['cell_line']=='PC9C1.1')] #Must be a pc9 cell line
t_pc9 = t_pc9[t_pc9['MCMC_converge']==1]     # Must converge 
t_pc9 = t_pc9[t_pc9['R2']>.5]    
#Remove selected combinations
t_pc9 = t_pc9[~np.logical_and(t_pc9['drug1_name']=='panobinostat',t_pc9['expt']=='HTS015_017_Combined.csv')]
t_pc9 = t_pc9[~np.logical_and(t_pc9['drug2_name']=='dasatinib',t_pc9['expt']=='HTS018_rates.csv')]
t_pc9 = t_pc9[~np.logical_and(t_pc9['drug2_name']=='trametinib',t_pc9['expt']=='HTS018_rates.csv')]
t_pc9 = t_pc9[~np.logical_and(t_pc9['drug2_name']=='bosutinib',t_pc9['expt']=='HTS018_rates.csv')]
#Must converge to model tier 4 or above
t_pc9 = t_pc9[t_pc9['model_level']>4]
#Enumerate through drug combination names
drug_combinations = t_pc9[['drug1_name','drug2_name']]      
#All drugs are in combination with osimertinib....Find the name of the other drug
drug_name=[]
for e,i in enumerate(drug_combinations.index):
    if drug_combinations['drug1_name'][i]!='osimertinib':
        drug_name.append(drug_combinations['drug1_name'][i])
    else:
        drug_name.append(drug_combinations['drug2_name'][i])                 
#Remove drugs with user noted error
bool_array = ~np.in1d(drug_name,['everolimus',#Darren noted error
'paclitaxel',#Darren noted error
'selumetinib',#Darren noted error
'buparlisib',#Darren noted error
'crizotinib',#Darren noted error
'alisertib',#Darren noted error
'erlotinib',#Darren noted error
'sunitinib',#Fluorescent
'cisplatin',#DMSO interaction
'vinorelbine',#duplicate
'pimobendan',#poorly classified mechanism
'5iodotubericidin',#poorly classified mechanism
'ml9',#poorly classified mechanism
'primaquinediphosphate'#poorly classified mechanism
])
t_pc9 = t_pc9[bool_array]
drug_name = list(compress(drug_name,bool_array))
drug_combinations = drug_combinations[bool_array]
#All drugs are in combination with osimertinib....Find the name of the other drug
drug_name=[]
for e,i in enumerate(drug_combinations.index):
    if drug_combinations['drug1_name'][i]!='osimertinib':
        drug_name.append(drug_combinations['drug1_name'][i])
    else:
        drug_name.append(drug_combinations['drug2_name'][i])  
#Write filtered results to csv...        
t_pc9.to_csv('../../Data/MasterResults_PC9_filtered.csv')

############################################################################        
#Read in table of associated mechanism for each drug 
############################################################################
mech = pd.read_csv('../../Data/Drug_MOA_Table_supp.csv',index_col=0)
mech_key = pd.read_csv('../../Data/Drug_Class_Key.csv',index_col=0)
#Which ones to label in the DSDs
to_lab = ['vindesine','vinorelbine','quisinostat',
          'm344','entinostat',
          'linsitinib','ceritinib']
#to_lab = drug_name
#Number of mechanism classes
num_classes=4
#Name of super classes
sup_class = ['Mitotic\nCheckpoint', 'Epigenetic\nRegulators','Receptors&\nChannels', 'Kinases']
mech = mech.loc[drug_name]
#Markers for each sub class
markers = ['D','s','o','^']
#Colors for each super class
col = ["blue", "orange", "red", "purple"]
cmaps = ["Blues","Oranges","Reds","Purples"]



################################################################################
###Plotting Data
#################################################################################
#Alpha values (X potentiates osimertinib)
x = [];x_sd=[]
for e,i in enumerate(t_pc9.index):
    if t_pc9['drug1_name'].loc[i]!='osimertinib':
        x.append(t_pc9['log_alpha1'].loc[i])
        x_sd.append(t_pc9['log_alpha1_std'].loc[i])
    else:
        x.append(t_pc9['log_alpha2'].loc[i])
        x_sd.append(t_pc9['log_alpha2_std'].loc[i])

x = np.array(x)
x_sd = np.array(x_sd)
#Beta values
y = np.array(t_pc9['beta_obs_norm'])
y_sd = np.array(t_pc9['beta_obs_norm_std'])
#Format axes
xmin,xmax = min(x)-.5,max(x)+.5
ymin,ymax = min(y)-.05,max(y)+.1

################################################################################
###Make initial density plot of the all the drugs
#################################################################################
fig = plt.figure(figsize=(1.2,1.2))
ax = [];
ax.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=3))
ax.append(plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=3))
ax.append(plt.subplot2grid((4,4),(1,3),rowspan=3,colspan=1))
plt.subplots_adjust(wspace=0,hspace=0)

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

ax[0].set_xlim((xmin,xmax))
ax[0].set_ylim((ymin,ymax))
ax[0].set_xlabel(r'$log(\alpha_1)$')
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
plt.savefig('Plots/Fig2_density_All.pdf')

################################################################################
###Make density plot for each super and sub classes
#################################################################################
#Density maps for super classes
for i in range(num_classes):
    fig = plt.figure(figsize=(1.2,1.2))
    ax = [];
    ax.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=3))
    ax.append(plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=3))
    ax.append(plt.subplot2grid((4,4),(1,3),rowspan=3,colspan=1))
    plt.subplots_adjust(wspace=0,hspace=0)
    a1 =  [x[drug_name.index(ind)] for ind in mech[(np.floor(mech['label'])==i+1)].index]   
    a2 =  [y[drug_name.index(ind)] for ind in mech[(np.floor(mech['label'])==i+1)].index]       
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
    ax[1].set_xticklabels([])
    ax[1].xaxis.set_label_position('top')
    ax[1].set_xlabel(('p-val:%.2f')%ks_2samp(a1,x)[1])
#    ax[1].set_xticklabels([str(int(j*100))+'%' for j in tiers])    
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
    ax[2].set_yticklabels([])
#    ax[2].set_yticklabels([str(int(j*100))+'%' for j in tiers])
    ax[2].plot(allD_by,allD_bx,'k',zorder=10)
    ax[2].yaxis.set_label_position('right')
    ax[2].set_ylabel(('p-val:%.2f')%ks_2samp(a2,y)[1])
    plt.savefig('Plots/Fig2_density_' + sup_class[i]+'.pdf')

#Density plots for sub-classes
for i in range(num_classes):
    for sub in mech['label'].loc[mech[(np.floor(mech['label'])==i+1)].index].unique():
        fig = plt.figure(figsize=(1.2,1.2))
        ax = [];
        ax.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=3))
        ax.append(plt.subplot2grid((4,4),(0,0),rowspan=1,colspan=3))
        ax.append(plt.subplot2grid((4,4),(1,3),rowspan=3,colspan=1))
        plt.subplots_adjust(wspace=0,hspace=0)
        a1 =  [x[drug_name.index(ind)] for ind in mech[mech['label']==sub].index]   
        a2 =  [y[drug_name.index(ind)] for ind in mech[mech['label']==sub].index] 
        if len(a1)>2:
            values = np.vstack([a1, a2])
            #Make a 2D gaussain
            kernel = st.gaussian_kde(values)
            f = np.reshape(kernel(positions).T, xx.shape)
            #intrp = interp2d(xx, yy, f, kind='linear', copy=True, bounds_error=False, fill_value=np.nan)
            tiers = [.01,.05,1]
            cfset = ax[0].contourf(xx, yy, f,tiers, cmap=cmaps[i])
            cset = ax[0].contour(xx, yy, f,tiers, colors='k')
            ax[0].clabel(cset, inline=1, fontsize=6)
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
            #ax[1].set_xticklabels([str(int(j*100))+'%' for j in tiers])
            ax[1].set_xticklabels([])
            ax[1].xaxis.set_label_position('top')
            ax[1].set_xlabel(('p-val:%.2f')%ks_2samp(a1,x)[1])
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
            ax[2].set_yticklabels([])
            ax[2].yaxis.set_label_position('right')
            ax[2].set_ylabel(('p-val:%.2f')%ks_2samp(a2,y)[1])
            ax[2].plot(allD_by,allD_bx,'k',zorder=10)
            plt.savefig('Plots/SupFig2_density_'+sup_class[i]+'_'+mech_key['subclass_to_write'].loc[sub]+'.pdf')


############################################################################
#Plot drug synergy diagrams for all, as well as each super and sub classes
############################################################################
#DSDs...
fig = plt.figure(figsize=(8.5,2))
ax = []
#num_classes+1 plots to include all subplot
for i in range(num_classes+1):
    ax.append(plt.subplot2grid((1,num_classes+1),(0,i)))
plt.subplots_adjust(wspace=0,hspace=.4)
#First figure is all drugs
plt.sca(ax[0])
plt.scatter(x,y,s=30,c='k',marker='o')
plt.errorbar(x,y,xerr=x_sd,yerr=y_sd,capsize=0,color='gray',lw=1,linestyle='')
plt.title('All')
ax[0].set_xlim((xmin,xmax))
ax[0].set_ylim((ymin,ymax))
y_lim = (ymin,ymax)
x_lim = (xmin,xmax)
plt.plot([0,0],y_lim,linestyle='--',color='k')
plt.plot(x_lim,[0,0],linestyle='--',color='k')
ax[0].tick_params(direction='in')
plt.xlabel(r'$log(\alpha_2)$ (X potentiates osi)')
plt.ylabel(r'$\bar{\beta_{obs}}$')

#Empty list of lists to hold the text labels 
text = [ [] for i in range(num_classes)]
for i in range(num_classes):
    plt.sca(ax[i+1])
    seen = []
    for ind in mech[(np.floor(mech['label'])==i+1)].index:
        if (mech['label'].loc[ind]!=0) & (np.in1d(ind,drug_name)[0]):
            mk = mech_key['subclass'].loc[mech['label'].loc[ind]]
            #If this type of marker has not been seen before, label point for 
            #legend and add plot title
            if ~np.in1d(mk,seen)[0]:
                seen.append(mk)
                plt.scatter(x[drug_name.index(ind)],y[drug_name.index(ind)],
                            s=30,color=col[i],marker=markers[int(np.round((mech['label'].loc[ind]-i-1)*10))-1],
                            zorder=10,edgecolor='k',linewidth=1,label=mk)
                plt.errorbar(x[drug_name.index(ind)],y[drug_name.index(ind)],xerr=x_sd[drug_name.index(ind)],yerr=y_sd[drug_name.index(ind)],capsize=0,color='gray',lw=1,linestyle='',zorder=100)

                plt.title(sup_class[i])
                #Add text label
                if np.in1d(ind,to_lab)[0]:                
                    text[i].append(plt.text(x[drug_name.index(ind)],y[drug_name.index(ind)],ind))
            #If the marker has been seen, dont label...
            else:
                plt.scatter(x[drug_name.index(ind)],y[drug_name.index(ind)],
                            s=30,color=col[i],marker=markers[int(np.round((mech['label'].loc[ind]-i-1)*10))-1],
                            zorder=10,edgecolor='k',linewidth=1)
                plt.errorbar(x[drug_name.index(ind)],y[drug_name.index(ind)],xerr=x_sd[drug_name.index(ind)],yerr=y_sd[drug_name.index(ind)],capsize=0,color='gray',lw=1,linestyle='',zorder=100)

                if np.in1d(ind,to_lab)[0]:                
                    text[i].append(plt.text(x[drug_name.index(ind)],y[drug_name.index(ind)],ind))
#Format all the axes
for e in range(num_classes):
    ax[e+1].set_xlim((xmin,xmax))
    ax[e+1].set_ylim((ymin,ymax))
    y_lim = (ymin,ymax)
    x_lim = (xmin,xmax)
    plt.sca(ax[e+1])
    plt.plot([0,0],y_lim,linestyle='--',color='k')
    plt.plot(x_lim,[0,0],linestyle='--',color='k')
    plt.sca(ax[e+1])
    #plt.legend(loc=9, bbox_to_anchor=(0.5, -0.03), ncol=1)
    plt.legend(loc='upper right',handletextpad=.1)
    ax[e+1].tick_params(direction='in')
    ax[e+1].set_xticklabels([])
    ax[e+1].set_yticklabels([])
    ax[e+1].xaxis.set_ticks_position('bottom')
    ax[e+1].yaxis.set_ticks_position('left')

#Adjust text function from 
#pip install adjustText package:  https://github.com/Phlya/adjustText
for i in range(num_classes):
    plt.sca(ax[i+1])
    adjust_text(text[i], arrowprops=dict(arrowstyle='->', color='red'))
plt.savefig('Fig2_DSDs.pdf',bbox_inches = 'tight',pad_inches = 0)







#Now plot the dsds for each of the sub classes
for i in range(num_classes):
    sub_classes = mech['label'].loc[mech[(np.floor(mech['label'])==i+1)].index].unique()
    text = [ [] for k in range(len(sub_classes)) ]
    fig = plt.figure(figsize=(1.7*len(sub_classes),1.5))
    ax = []
    for k in range(len(sub_classes)):
        ax.append(plt.subplot2grid((1,len(sub_classes)),(0,k)))
    plt.subplots_adjust(wspace=0,hspace=.4)
    for e,sub in enumerate(sub_classes):
        plt.sca(ax[e])
        seen = []
        for ind in mech[mech['label']==sub].index:
            mk = mech_key['subclass'].loc[mech['label'].loc[ind]]
            plt.scatter(x[drug_name.index(ind)],y[drug_name.index(ind)],
                        s=30,c=col[i],marker=markers[int(np.round((mech['label'].loc[ind]-i-1)*10))-1],
                        zorder=10,edgecolor='k',linewidth=1)
            plt.title(mech_key['subclass'].loc[sub])
            if np.in1d(ind,to_lab)[0]:                
                text[e].append(plt.text(x[drug_name.index(ind)],y[drug_name.index(ind)],ind))

        
    for e in range(len(sub_classes)):
        ax[e].set_xlim((xmin,xmax))
        ax[e].set_ylim((ymin,ymax))
        y_lim = (ymin,ymax)
        x_lim = (xmin,xmax)
        plt.sca(ax[e])
        plt.plot([0,0],y_lim,linestyle='--',color='k')
        plt.plot(x_lim,[0,0],linestyle='--',color='k')
        plt.sca(ax[e])
        #plt.legend(loc=9, bbox_to_anchor=(0.5, -0.03), ncol=1)
        plt.legend(loc='upper right')
        ax[e].tick_params(direction='in')
        if e!=0:
            ax[e].set_xticklabels([])
            ax[e].set_yticklabels([])
    plt.sca(ax[0])    
    plt.xlabel(r'$log(\alpha_2)$ (X potentiates osi)')
    plt.ylabel(r'$\bar{\beta_{obs}}$')
    for e in range(len(sub_classes)):
        plt.sca(ax[e])
        adjust_text(text[e], arrowprops=dict(arrowstyle='->', color='red'))
    plt.savefig('Plots/SupFig2_DSDs_'+sup_class[i]+'_subClasses.pdf')


#########################################################################################
#Now make figure labeling each drugs super mechanism 
#########################################################################################
fig = plt.figure(figsize=(7.5,2))
ax1 = plt.subplot(111)
ax1.set_ylim((0,num_classes))
ax1.set_xlim((0,len(drug_name)+1))
#Set up ticks and grid
ax1.set_xticks(np.linspace(0.5,len(drug_name)-.5,len(drug_name)))
ax1.set_yticks(np.linspace(0.5,num_classes-.5,num_classes))
ax1.set_xticks(np.linspace(0,len(drug_name),len(drug_name)+1),minor=True)
ax1.set_yticks(np.linspace(0,num_classes-1,num_classes),minor=True)
ax1.tick_params(which='minor',length=0)
plt.grid(axis='both',color='k',linewidth=1,which='minor')
plt.grid(axis='both',color='k',linewidth=1,which='minor')
ax1.xaxis.tick_top()
#Sort drug names alphabetically
idd = np.argsort(drug_name)
#drug_name[drug_name.index('beclomethasone')]='beclomethasonedipropionate'
#Plot a circle
for e,i in enumerate(t_pc9.index):
    if t_pc9['drug1_name'].loc[i]=='osimertinib':
        il = 2
    else:
        il = 1
    idx = float(np.floor(mech['label'].loc[t_pc9['drug'+str(il)+'_name'].loc[i]]))
    idx2 = drug_name.index(t_pc9['drug'+str(il)+'_name'].loc[i])
    if idx!=0:
        ax1.add_artist(plt.Circle((np.where(idd==idx2)[0][0]+.5,idx-.5),radius=.5,color=col[int(idx-1)]))        
drug_name[drug_name.index('beclomethasonedipropionate')]='beclomethasone'
ax1.set_xticklabels([drug_name[i] for i in idd],fontsize=6)
plt.setp(ax1.xaxis.get_majorticklabels(), rotation=-45 ) 
plt.xticks(ha='right')
ax1.set_yticklabels(sup_class,fontsize=8,ha='right')
plt.tight_layout()
plt.savefig('Plots/Fig2_drugClasses.pdf')

###################################################################################
#####Make boxplot of E1,E2,E3 for each super class and sub class
###################################################################################
#fig = plt.figure(figsize=(8.5,1.5))
#ax = []
#for i in range(num_classes+1):
#    ax.append(plt.subplot2grid((1,num_classes+1),(0,i)))
#plt.subplots_adjust(wspace=0,hspace=.4)
#text = [ [] for i in range(num_classes)]
#mech = mech.loc[drug_name]
#e1 = [] #Drug X effect
#for i in t_pc9.index:
#    if t_pc9['drug1_name'].loc[i]!='osimertinib':
#        e1.append(t_pc9['E1_obs'].loc[i])
#    else:
#        e1.append(t_pc9['E2_obs'].loc[i])
#e2 = [] #Osimertinib max effect
#for i in t_pc9.index:
#    if t_pc9['drug1_name'].loc[i]!='osimertinib':
#        e2.append(t_pc9['E2_obs'].loc[i])
#    else:
#        e2.append(t_pc9['E1_obs'].loc[i])
#e1 = np.array(e1)
#e2 = np.array(e2)                
#e3 = t_pc9['E3_obs'].values #Combined max effect
#
##Make plot for all drugs
#plt.sca(ax[0])
#x1 = 0+(np.random.rand(len(e1))-.5)/5
#x2 = 1+(np.random.rand(len(e2))-.5)/5
#x3 = 2+(np.random.rand(len(e3))-.5)/5
#plt.scatter(x1,e2,s=30,c='k',marker='o',zorder=10,edgecolor='k',linewidth=1)
#plt.scatter(x2,e1,s=30,c='k',marker='o',zorder=10,edgecolor='k',linewidth=1)
#plt.scatter(x3,e3,s=30,c='k',marker='o',zorder=10,edgecolor='k',linewidth=1)
#plt.title('All')            
#for j in range(len(e1)):
#    plt.plot((x2[j],x3[j]),(e1[j],e3[j]),linewidth=.5,color='k')
#bp = plt.boxplot([e2,e1,e3],positions=[0,1,2],showfliers=False,zorder=100)
#ax[0].set_xlim([-.25,2.25])            
#ax[0].set_ylim((min(e3),t_pc9['E0'].max()))
#ax[0].set_xticks([0,1,2])
#ax[0].set_ylabel('DIP s-1')
#ax[0].set_xticklabels(['Max Osi\n(E1)','Max X\n(E2)','Max Osi+X\n(E3)'])
#plt.plot([-.25,2.25],[0,0],linestyle='--',color='k')        
#plt.plot([-.25,2.25],[np.mean(e2),np.mean(e2)],linestyle='-',color='r',label='Max Osi',zorder=1000)
#plt.legend(loc='lower left')
#
#
##make plot for each super class only plotting max drug X (E2) and max drug X+max drug osi (E3)
#for i in range(num_classes):
#    plt.sca(ax[i+1])
#    seen = []
#    a1 = [e1[drug_name.index(j)] for j in mech[(np.floor(mech['label'])==i+1)].index]
#    a2 = [e2[drug_name.index(j)] for j in mech[(np.floor(mech['label'])==i+1)].index]
#    a3 = [e3[drug_name.index(j)] for j in mech[(np.floor(mech['label'])==i+1)].index]
#    x1 = 0+(np.random.rand(len(a1))-.5)/5
#    x2 = 1+(np.random.rand(len(a2))-.5)/5
#    x3 = 2+(np.random.rand(len(a3))-.5)/5
#    text = []
#    for en,ind in enumerate(mech[(np.floor(mech['label'])==i+1)].index):
#        plt.scatter(x2[en],e1[drug_name.index(ind)],
#                    s=30,c=col[i],marker=markers[int(np.round((mech['label'].loc[ind]-i-1)*10))-1],
#                    zorder=10,edgecolor='k',linewidth=1)
#        plt.scatter(x3[en],e3[drug_name.index(ind)],
#                    s=30,c=col[i],marker=markers[int(np.round((mech['label'].loc[ind]-i-1)*10))-1],
#                    zorder=10,edgecolor='k',linewidth=1)
#        if np.in1d(ind,to_lab)[0]:                
#            text.append(plt.text(x3[en],e3[drug_name.index(ind)],ind))
#    plt.title(sup_class[i])
#    for j in range(len(a2)):
#        plt.plot((x2[j],x3[j]),(a1[j],a3[j]),linewidth=.5,color='k')
#    
#    bp = plt.boxplot([a1,a3],positions=[1,2],showfliers=False,zorder=10000,boxprops = {'linewidth':2},whiskerprops={'linewidth':2},medianprops = {'linewidth':2})
#    ax[i+1].set_xlim([.75,2.25])            
#    ax[i+1].set_ylim((min(e3),t_pc9['E0'].max()))
#    ax[i+1].set_xticks([1,2])
#    plt.plot([.75,2.25],[0,0],linestyle='--',color='k')
#    plt.plot([.75,2.25],[np.mean(e2),np.mean(e2)],linestyle='-',color='r',zorder=100000)
#    if i==0:
#        ax[i+1].set_xticklabels(['E2','E3'])
#        ax[i+1].set_yticks([])
#    else:
#        ax[i+1].set_xticklabels([])
#        ax[i+1].set_yticks([])
#    adjust_text(text, arrowprops=dict(arrowstyle='->', color='red'))
#
#plt.savefig('Plots/Fig2_effect.pdf')
#
#
#
##Now for all the sub-classes....
#for i in range(num_classes):
#    sub_classes = mech['label'].loc[mech[(np.floor(mech['label'])==i+1)].index].unique()
#    text = [ [] for k in range(len(sub_classes)) ]
#    fig = plt.figure(figsize=(1.5*len(sub_classes),1.5))
#    ax = []
#    for k in range(len(sub_classes)):
#        ax.append(plt.subplot2grid((1,len(sub_classes)),(0,k)))
#    text = [ [] for k in range(len(sub_classes)) ]
#    plt.subplots_adjust(wspace=0,hspace=.4)
#    for e,sub in enumerate(sub_classes):
#        plt.sca(ax[e])
#        seen = []
#        a1 = [e1[drug_name.index(j)] for j in mech[mech['label']==sub].index]
#        a2 = [e2[drug_name.index(j)] for j in mech[mech['label']==sub].index]
#        a3 = [e3[drug_name.index(j)] for j in mech[mech['label']==sub].index]
#        x1 = 0+(np.random.rand(len(a1))-.5)/5
#        x2 = 1+(np.random.rand(len(a2))-.5)/5
#        x3 = 2+(np.random.rand(len(a3))-.5)/5
#        for en,ind in enumerate(mech[mech['label']==sub].index):
#            plt.scatter(x2[en],e1[drug_name.index(ind)],
#                        s=30,c=col[i],marker=markers[int(np.round((mech['label'].loc[ind]-i-1)*10))-1],
#                        zorder=10,edgecolor='k',linewidth=1)
#            plt.scatter(x3[en],e3[drug_name.index(ind)],
#                        s=30,c=col[i],marker=markers[int(np.round((mech['label'].loc[ind]-i-1)*10))-1],
#                        zorder=10,edgecolor='k',linewidth=1)
#            if np.in1d(ind,to_lab)[0]:                
#                text[e].append(plt.text(x3[en],e3[drug_name.index(ind)],ind))
#       
#        for j in range(len(a2)):
#            plt.plot((x2[j],x3[j]),(a1[j],a3[j]),linewidth=.5,color='k')
#        
#        bp = plt.boxplot([a1,a3],positions=[1,2],showfliers=False,zorder=10000,boxprops = {'linewidth':2},whiskerprops={'linewidth':2},medianprops = {'linewidth':2})
#        ax[e].set_xlim([.75,2.25])            
#        ax[e].set_ylim((min(e3),t_pc9['E0'].max()))
#        ax[e].set_xticks([1,2])
#        plt.plot([.75,2.25],[0,0],linestyle='--',color='k')
#        plt.plot([.75,2.25],[np.mean(e2),np.mean(e2)],linestyle='-',color='r',zorder=10001)
#        if e==0:
#            ax[e].set_xticklabels(['E2','E3'])
#            ax[e].set_yticks([])
#        else:
#            ax[e].set_xticklabels([])
#            ax[e].set_yticks([])
#        adjust_text(text[e], arrowprops=dict(arrowstyle='->', color='red'))
#        plt.title(mech_key['subclass'].loc[sub])
#        plt.savefig('Plots/SupFig2_effect_'+sup_class[i]+'.pdf')


############################################################################################################################################
####Plot static surfaces
############################################################################################################################################
os.chdir('../Functions/')
import surfPlots as sP
os.chdir('../Fig2/')
for drug, drug_id in zip(["m344","ceritinib"], ["drug1","drug2"]):
    drug = drug.lower()

    hts_data_18 = pd.read_csv("../../Data/HTS018_rates.csv")
    drug1_name = [i.lower() for i in hts_data_18['drug1']]
    drug2_name = [i.lower() for i in hts_data_18['drug2']]
    hts_data_18['drug1'] = drug1_name
    hts_data_18['drug2'] = drug2_name

    hts_data_15 = pd.read_csv("../../Data/HTS015_017_Combined.csv")
    drug1_name = [i.lower() for i in hts_data_15['drug1']]
    drug2_name = [i.lower() for i in hts_data_15['drug2']]
    hts_data_15['drug1'] = drug1_name
    hts_data_15['drug2'] = drug2_name

    hts_data_15 = hts_data_15.loc[((hts_data_15['drug1']==drug) | (hts_data_15['drug1.conc']==0)) & (hts_data_15['cell.line']=="PC9c1")]
    hts_data_18 = hts_data_18.loc[((hts_data_18['drug2']==drug) | (hts_data_18['drug2.conc']==0)) & (hts_data_18['cell.line']=="PC9c1")]

    if drug_id == "drug1": hts_data = hts_data_15
    else: hts_data = hts_data_18

    mcmc_data = pd.read_csv("../../Data/MasterResults_PC9_filtered.csv")
    mcmc_data = mcmc_data.loc[mcmc_data['%s_name'%drug_id]==drug]

    sP.plot_surface(mcmc_data, hts_data, drug_id, drug, fname="%s.pdf"%drug, zlim = (-0.03, 0.032))
    sP.plot_slices(mcmc_data, hts_data, drug_id, drug, fname="%s_slices.pdf"%drug)


##################################################
###### Make interactive HTML plot in html file for all combinations
##################################################
df = pd.read_csv('../../Data/MasterResults_PC9_filtered.csv')
for idx in df.index:
    ##Now plot the dose response slices...
    T = df.loc[idx].to_dict()
        
    #Should the dose response surfaces be plotted?
    to_plot = 1#Should the results be plotted?
    expt = '../../Data/' + T['expt']#Experiment file
    target_cell_line = T['cell_line']#The cell line to consider
    
    #Read in the data into a pandas data frame
    data = pd.read_table(expt, delimiter=',')
    data['cell.line'] = data['cell.line'].str.upper()
    data['drug1'] = data['drug1'].str.lower()
    data['drug2'] = data['drug2'].str.lower()
    #Subset by target cell line
    data = data[data['cell.line']==T['cell_line']]
    #data = data.reset_index()
    
    drug1_name = T['drug1_name']#Drug1 name to replace in array gen function
    drug2_name = T['drug2_name']#Drug2 name to replace in ArrayGen Function
    
    tmp_data = data[
                    ((data['drug1']=='control') & (data['drug2']=='control'))   | 
                    ((data['drug1']==drug1_name) & (data['drug2.conc']==0))    | 
                    ((data['drug2']==drug2_name) & (data['drug1.conc']==0))     |
                    ((data['drug2']==drug2_name) & (data['drug1']==drug1_name))
                    ]
    
    #tmp_data = tmp_data.reset_index() #Reset the data frame index
    d1 = tmp_data['drug1.conc'].values #Vector of drug 1 concentrations
    d2 = tmp_data['drug2.conc'].values #Vector of drug2 concentrations
    dip = tmp_data['rate'].values #Vector of the measured dip rate
    dip_95ci = tmp_data['rate.95ci'].values #Vector of the 95% confidence interval on the DIP fit.
    dip_sd = dip_95ci/(2*1.96)
        
    #Remove nan values
    d1      = d1[~np.isnan(dip)]
    d2      = d2[~np.isnan(dip)]
    dip_sd  = dip_sd[~np.isnan(dip)]
    dip     = dip[~np.isnan(dip)]
                  
    
    #Force all the drug names to lower case...
    drug1_name = drug1_name.lower()
    drug2_name = drug2_name.lower()
    #Force all the cell_line names to be upper case
    target_cell_line = target_cell_line.upper()
    
    
    mod_lev = T['model_level']           
    #If model level was less than 5, then it is some approximation of the Detail balance case.
    if mod_lev < 5:     
        popt3 = [ T['E0'],
                  T['E0_std'],
                  T['E1'],
                  T['E1_std'],
                  T['E2'],
                  T['E2_std'],
                  T['E3'],
                  T['E3_std'],
                  10**T['log_C1'],
                  np.log(10)*T['log_C1_std']*10**T['log_C1'],
                  10**T['log_C2'],
                  np.log(10)*T['log_C2_std']*10**T['log_C2'],
                  T['h1'],
                  T['h1_std'],
                  T['h2'],
                  T['h2_std'],
                  10**T['log_alpha1'],
                  np.log(10)*T['log_alpha1_std']*10**T['log_alpha1']]
        sP.DosePlots_PLY(d1,d2,dip,dip_sd,drug1_name,drug2_name,popt3[0:17:2],'DB',target_cell_line,expt,'DIP rate (h-1)')
    elif mod_lev == 5: #It was best fit by a model which did not obey detail balance
        popt3 = [ T['E0'],
                  T['E0_std'],
                  T['E1'],
                  T['E1_std'],
                  T['E2'],
                  T['E2_std'],
                  T['E3'],
                  T['E3_std'],
                  T['r1'],
                  T['r1_std'],
                  T['r2'],
                  T['r2_std'],
                  10**T['log_C1'],
                  np.log(10)*T['log_C1_std']*10**T['log_C1'],
                  10**T['log_C2'],
                  np.log(10)*T['log_C2_std']*10**T['log_C2'],
                  T['h1'],
                  T['h1_std'],
                  T['h2'],
                  T['h2_std'],
                  10**T['log_alpha1'],
                  np.log(10)*T['log_alpha1_std']*10**T['log_alpha1'],
                  10**T['log_alpha2'],
                  np.log(10)*T['log_alpha2_std']*10**T['log_alpha2']]
        sP.DosePlots_PLY(d1,d2,dip,dip_sd,drug1_name,drug2_name,popt3[0:23:2],'NDB',target_cell_line,expt,'DIP rate (h-1)')        
        
    else:
        popt3 = [ T['E0'],
                  T['E0_std'],
                  T['E1'],
                  T['E1_std'],
                  T['E2'],
                  T['E2_std'],
                  T['E3'],
                  T['E3_std'],
                  T['r1'],
                  T['r1_std'],
                  T['r2'],
                  T['r2_std'],
                  10**T['log_C1'],
                  np.log(10)*T['log_C1_std']*10**T['log_C1'],
                  10**T['log_C2'],
                  np.log(10)*T['log_C2_std']*10**T['log_C2'],
                  T['h1'],
                  T['h1_std'],
                  T['h2'],
                  T['h2_std'],
                  10**T['log_alpha1'],
                  np.log(10)*T['log_alpha1_std']*10**T['log_alpha1'],
                  10**T['log_alpha2'],
                  np.log(10)*T['log_alpha2_std']*10**T['log_alpha2'],
                  10**T['log_gamma1'],
                  np.log(10)*T['log_gamma1_std']*10**T['log_gamma1'],
                  10**T['log_gamma2'],
                  np.log(10)*T['log_gamma2_std']*10**T['log_gamma2']]
        sP.DosePlots_PLY(d1,d2,dip,dip_sd,drug1_name,drug2_name,popt3[0:27:2],'NDB_hill',target_cell_line,expt,'DIP rate (h-1)')          


