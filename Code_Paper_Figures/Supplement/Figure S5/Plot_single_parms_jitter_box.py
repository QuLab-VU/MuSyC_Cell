#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:01:49 2017

@author: xnmeyer
"""

import pandas as pd
import numpy as np
from itertools import compress
import matplotlib.pyplot as plt
from matplotlib import rc
font = {'family' : 'normal',
        'weight':'normal',
        'size'   : 8}
axes = {'linewidth': 3}
rc('font', **font)
rc('axes',**axes)
from scipy.integrate import quad
from matplotlib.artist import setp
t_pc9 = pd.read_csv('../../../Data/MasterResults_PC9_filtered.csv')
drug_combinations = t_pc9[['drug1_name','drug2_name']]      

#All drugs are in combination with osimertinib....Find the name of the other drug
drug_name=[]
for e,i in enumerate(drug_combinations.index):
    if drug_combinations['drug1_name'][i]!='osimertinib':
        drug_name.append(drug_combinations['drug1_name'][i])
    else:
        drug_name.append(drug_combinations['drug2_name'][i])   
        
#Calculate all the different typ
def Edrug1D(d,h,E0,Em,C):
    return Em + (E0-Em) / (1 + (d/C)**h)

def integrand(x,a,b,c,d):
    return d+(a-d)/(1+(x/c)**b)
    
E_obs = np.zeros((len(t_pc9),2))
E_inf = np.zeros((len(t_pc9),2))
E_50 = np.zeros((len(t_pc9),2))
AUC = np.zeros((len(t_pc9),2))
h = np.zeros((len(t_pc9),2))
drug_name = []  
for e,i in enumerate(t_pc9.index):
        E_obs[e,0]=Edrug1D(float(t_pc9['max_conc_d1'][i]),float(t_pc9['h1'][i]),float(t_pc9['E0'][i]),float(t_pc9['E1'][i]),10**float(t_pc9['log_C1'][i]))
        E_obs[e,1]=Edrug1D(float(t_pc9['max_conc_d2'][i]),float(t_pc9['h2'][i]),float(t_pc9['E0'][i]),float(t_pc9['E2'][i]),10**float(t_pc9['log_C2'][i]))
        E_inf[e,0] = float(t_pc9['E1'][i])
        E_inf[e,1] = float(t_pc9['E2'][i])
        E_50[e,0] = 10**float(t_pc9['log_C1'][i])
        E_50[e,1] = 10**float(t_pc9['log_C2'][i])
        
        h[e,0] = float(t_pc9['h1'][i])
        h[e,1] = float(t_pc9['h2'][i])
        AUC[e,0] = quad(integrand,  float(t_pc9['max_conc_d1'][i]),
                                    float(t_pc9['min_conc_d1'][i]), 
                                    args=(float(t_pc9['E0'][i]),
                                          float(t_pc9['h1'][i]),
                                          10**float(t_pc9['log_C1'][i]),
                                          float(t_pc9['E1'][i])))[0]
        AUC[e,1] = quad(integrand,  float(t_pc9['max_conc_d2'][i]),
                                    float(t_pc9['min_conc_d2'][i]), 
                                    args=(float(t_pc9['E0'][i]),
                                          float(t_pc9['h2'][i]),
                                          10**float(t_pc9['log_C2'][i]),
                                          float(t_pc9['E2'][i])))[0]
        drug_name.append(t_pc9['drug1_name'][i])
        drug_name.append(t_pc9['drug2_name'][i])


E_obs = np.concatenate([E_obs[:,0],E_obs[:,1]])
h = np.concatenate([h[:,0],h[:,1]])
E_inf = np.concatenate([E_inf[:,0],E_inf[:,1]])
E_50 = np.concatenate([E_50[:,0],E_50[:,1]])
AUC = np.concatenate([AUC[:,0],AUC[:,1]])

df = pd.DataFrame({'$E_{max}obs$':E_obs,'drug':drug_name,'h':h,'$log(E_{50})$':np.log10(E_50)})
df = df.set_index('drug')
df = df.groupby(by='drug').mean()
drug_name = list(df.index)
to_label = ['osimertinib','vindesine','ceritinib','m344']

jitter_width = .4
fig = plt.figure(figsize=(3,4))


ax=[]
colors = ['r','g','m','y','c']
for e,parm in enumerate(list(df)):
    ax.append(plt.subplot(3,1,e+1))
    ax[e].spines['right'].set_color('none')
    ax[e].spines['top'].set_color('none')
    ax[e].spines['left'].set_color('none')
    pt = plt.boxplot(np.array(df[parm]),vert=False,widths=.5,showfliers=False)
    for i in pt.items():
        setp(pt[i[0]],linewidth=3)
    yl = plt.ylim()
    y = np.random.uniform(yl[0]+jitter_width*(yl[1]-yl[0]),yl[1]-jitter_width*(yl[1]-yl[0]),len(df[parm]))
    sc = [plt.scatter(np.array(df[parm]),y,edgecolors='k',linewidth=1,s=20)]
    #for tick in ax1.get_xticklabels():
    #    tick.set_rotation(45)
    plt.ylabel(parm)
    ax[e].set_yticks([])
    for label in to_label:
        xt = np.array(df[parm].loc[label])
        yt = y[drug_name==label]
        sc.append(plt.scatter(xt,yt,c=colors[to_label.index(label)],s=40,edgecolor='k',linewidth=1,label=label))
#        plt.annotate(
#                    label,
#                    xy=(xt, yt), xytext=(-10, 10),
#                    textcoords='offset points', ha='left', va='bottom')
    ax[e].patch.set_facecolor('None')
    ax[e].locator_params(nbins=4, axis='x')
    for e,pid in enumerate(sc):
        pid.set_zorder(20+e)

plt.sca(ax[2])
plt.legend(loc=9, bbox_to_anchor=(0.5, -0.15), ncol=len(to_label))
plt.xticks((1,5,10))
plt.tight_layout()
plt.savefig('SingleParms.pdf')



fig = plt.figure(figsize=(2.5,4))
ax = []
for e in range(3):
    ax.append(plt.subplot(3,1,e+1))
    ax[e].spines['right'].set_color('none')
    ax[e].spines['top'].set_color('none')
    ax[e].spines['bottom'].set_color('k')
    ax[e].spines['left'].set_color('k')
    plt.ylabel(df.columns[e])
    ax[e].set_facecolor('w')
    
x1=[];x2=[];x3=[];y1=[];y2=[];y3=[];
    
for e,i in enumerate(t_pc9.index):
    if t_pc9['drug1_name'].loc[i]=='osimertinib':
        idt = 2
    else:
        idt = 1
    plt.sca(ax[0])
    plt.scatter(t_pc9['beta_obs_norm'].loc[i],t_pc9['E'+str(idt)+'_obs'].loc[i],s=30,c='k')
    x1.append(t_pc9['beta_obs_norm'].loc[i])
    y1.append(t_pc9['E'+str(idt)+'_obs'].loc[i])
    plt.sca(ax[1])    
    plt.scatter(t_pc9['log_alpha'+str(idt)].loc[i],t_pc9['log_C'+str(idt)].loc[i],s=30,c='k')
    x2.append(t_pc9['log_alpha'+str(idt)].loc[i])
    y2.append(t_pc9['log_C'+str(idt)].loc[i])
    plt.sca(ax[2])
    plt.scatter(t_pc9['log_alpha'+str(idt)].loc[i],t_pc9['h'+str(idt)].loc[i],s=30,c='k')
    x3.append(t_pc9['log_alpha'+str(idt)].loc[i])
    y3.append(t_pc9['h'+str(idt)].loc[i])
    
    
plt.sca(ax[0])
plt.xlabel(r'$\beta_{obs}$')
plt.text(.75,.025,'Corr: ' + str("%.2f" % np.corrcoef(x1,y1)[0,1]),fontsize=8,bbox=dict(facecolor='none', edgecolor='k'))
plt.sca(ax[1])
plt.xlabel(r'$log(\alpha_2)$(X potentiates osi)')
plt.text(-4,-10,'Corr: ' + str("%.2f" % np.corrcoef(x2,y2)[0,1]),fontsize=8,bbox=dict(facecolor='none', edgecolor='k'))
plt.sca(ax[2])
plt.xlabel(r'$log(\alpha_2)$(X potentiates osi)')
plt.text(5,7,'Corr: ' + str("%.2f" % np.corrcoef(x3,y3)[0,1]),fontsize=8,bbox=dict(facecolor='none', edgecolor='k'))
plt.tight_layout()
plt.savefig('SingleParms_synCorr.pdf')
