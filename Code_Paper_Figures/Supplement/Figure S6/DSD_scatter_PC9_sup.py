# -*- coding: utf-8 -*-
"""
Christian T Meyer
email christian.t.meyer@vanderbilt.edu
Code for creating figure S3 in paper 'Quantifying Drug Synergy with respect to potency and efficacy'
Commented and up to date 2-12-18
"""
import pandas as pd
import numpy as np
from itertools import compress
import matplotlib.pyplot as plt
from matplotlib import rc
from adjustText import adjust_text
font = {'family' : 'normal',
        'weight':'normal',
        'size'   : 8}
axes = {'linewidth': 2}
rc('font', **font)
rc('axes',**axes)

t_pc9 = pd.read_csv('../../../Data/MasterResults_PC9_filtered.csv')
drug_combinations = t_pc9[['drug1_name','drug2_name']]      

#All drugs are in combination with osimertinib....Find the name of the other drug
drug_name=[]
for e,i in enumerate(drug_combinations.index):
    if drug_combinations['drug1_name'][i]!='osimertinib':
        drug_name.append(drug_combinations['drug1_name'][i])
    else:
        drug_name.append(drug_combinations['drug2_name'][i])   
#All drugs are in combination with osimertinib....Find the name of the other drug
drug_name=[]
for e,i in enumerate(drug_combinations.index):
    if drug_combinations['drug1_name'][i]!='osimertinib':
        drug_name.append(drug_combinations['drug1_name'][i])
    else:
        drug_name.append(drug_combinations['drug2_name'][i])   
        
#Read in table of associated mechanism for each drug 
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


x = []
for e,i in enumerate(t_pc9.index):
    if t_pc9['drug1_name'].loc[i]!='osimertinib':
        x.append(t_pc9['log_alpha2'].loc[i])
    else:
        x.append(t_pc9['log_alpha1'].loc[i])
x = np.array(x)
#Beta values
y = np.array(t_pc9['beta_obs_norm'])
#Format axes
xmin,xmax = min(x)-.5,max(x)+.5
ymin,ymax = min(y)-.05,max(y)+.1

#Now make figure with alpha alpha plot supplement
fig = plt.figure(figsize=(7.5,2))
ax = []
for i in range(num_classes):
    ax.append(plt.subplot2grid((1,num_classes),(0,i)))
plt.subplots_adjust(wspace=0,hspace=.4)

#Empty list of lists to hold the text labels 
text = [ [] for i in range(num_classes)]
for i in range(num_classes):
    plt.sca(ax[i])
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
    ax[e].set_xlim((xmin,xmax))
    ax[e].set_ylim((ymin,ymax))
    y_lim = (ymin,ymax)
    x_lim = (xmin,xmax)
    plt.sca(ax[e])
    plt.plot([0,0],y_lim,linestyle='--',color='k')
    plt.plot(x_lim,[0,0],linestyle='--',color='k')
    plt.sca(ax[e])
    #plt.legend(loc=9, bbox_to_anchor=(0.5, -0.03), ncol=1)
    plt.legend(loc='lower right')
    ax[e].tick_params(direction='in')
    ax[e].set_xticklabels([])
    ax[e].set_yticklabels([])
    
#Adjust text function from 
#pip install adjustText package:  https://github.com/Phlya/adjustText
for i in range(num_classes):
    plt.sca(ax[i])
    adjust_text(text[i], arrowprops=dict(arrowstyle='->', color='red'))
plt.sca(ax[0])    
plt.xlabel(r'$log(\alpha_1)$(X Potentiates Osi)')
plt.ylabel(r'$log(\alpha_2)$(Osi Potentiates X')
plt.savefig('DSDs_Supp_PC9.pdf')
