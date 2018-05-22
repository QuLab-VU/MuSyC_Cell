# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 15:42:58 2018

@author: xnmeyer
"""

#Import libraries and configure default font/axes parameters
import pandas as pd
import numpy as np
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
t_pc9 = pd.read_csv('../../../Data/MasterResults_PC9_filtered.csv')
#Enumerate through drug combination names
drug_combinations = t_pc9[['drug1_name','drug2_name']]      
#All drugs are in combination with osimertinib....Find the name of the other drug
drug_name=[]
for e,i in enumerate(drug_combinations.index):
    if drug_combinations['drug1_name'][i]!='osimertinib':
        drug_name.append(drug_combinations['drug1_name'][i])
    else:
        drug_name.append(drug_combinations['drug2_name'][i])                 
############################################################################        
#Read in table of associated mechanism for each drug 
############################################################################
mech = pd.read_csv('../../../Data/Drugs_MOA_3-22-18.csv',index_col=0)
mech_key = pd.read_csv('../../../Data/Drug_Class_Key.csv',index_col=0)
#Which ones to label in the DSDs
to_lab = ['vindesine','vinorelbine',
          'm344',
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
col = ['tab:blue', 'tab:orange', 'tab:red', 'tab:purple']
cmaps = ['Blues','Oranges','Reds','Purples']



################################################################################
###Make initial density plot of the all the drugs
#################################################################################
#Alpha values (X potentiates osimertinib)
x = []
for e,i in enumerate(t_pc9.index):
    if t_pc9['drug1_name'].loc[i]!='osimertinib':
        x.append(t_pc9['log_alpha2'].loc[i])
    else:
        x.append(t_pc9['log_alpha1'].loc[i])
x = np.array(x)
#Beta values
y = np.array(t_pc9['beta_obs_norm'])
xmin,xmax = min(x)-.5,max(x)+.5
ymin,ymax = min(y)-.05,max(y)+.1

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
                            s=30,c=col[i],marker=markers[int(np.round((mech['label'].loc[ind]-i-1)*10))-1],
                            zorder=10,edgecolor='k',linewidth=1,label=mk)
                plt.title(sup_class[i])
                #Add text label
                if np.in1d(ind,to_lab)[0]:                
                    text[i].append(plt.text(x[drug_name.index(ind)],y[drug_name.index(ind)],ind))
            #If the marker has been seen, dont label...
            else:
                plt.scatter(x[drug_name.index(ind)],y[drug_name.index(ind)],
                            s=30,c=col[i],marker=markers[int(np.round((mech['label'].loc[ind]-i-1)*10))-1],
                            zorder=10,edgecolor='k',linewidth=1)
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
    plt.legend(loc='upper right')
    ax[e+1].tick_params(direction='in')
    ax[e+1].set_xticklabels([])
    ax[e+1].set_yticklabels([])
    
#Adjust text function from 
#pip install adjustText package:  https://github.com/Phlya/adjustText
for i in range(num_classes):
    plt.sca(ax[i+1])
    adjust_text(text[i], arrowprops=dict(arrowstyle='->', color='red'))
plt.savefig('Plots/SuppFig_DSDs.pdf')


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
