#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 08:28:38 2018
Code for generating figure relating to fitting algorithm
"""
#Import packages
import h5py
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import rc
#Adjust text from  https://github.com/Phlya/adjustText
#Can be installed with pip
from adjustText import adjust_text
import scipy.stats as st
font = {'family' : 'arial',
        'weight':'normal',
        'size'   : 8}
axes = {'linewidth': 2}
rc('font', **font)
rc('axes',**axes)

#Data generated in geneterate_test_data.py.  
#After fitting move to top directory
dirs = glob.glob('data*/*/')
#Compile all the results
T = pd.DataFrame([])
alpha_vals = [-2,-1,0,1,2]
beta_vals = [-.5,-.25,0,.25,.5]
for d in dirs:
    try:
        tmp = pd.read_csv(d+'MasterResults.csv')
        if 'DIC_0' in list(tmp) and tmp['fit_algorithm'].unique()[0]=='pso_mcmc':
            tmp['fit_algorithm']='pso_multi_mcmc'
        alpha_indx = [int(i[2][5:]) for i in tmp['cell_line'].str.split('_')]
        beta_indx =  [int(i[3][4:]) for i in tmp['cell_line'].str.split('_')]
        beta = [beta_vals[i] for i in beta_indx]
        log_alpha = [alpha_vals[i] for i in alpha_indx]
        tmp['beta_true'] = np.array(beta)
        tmp['E0_true'] = .3
        tmp['E1_true']=0.
        tmp['E2_true']=0.
        tmp['E3_true']=-np.array(beta)*.3 
        tmp['h1_true']=1.
        tmp['h2_true']=1.
        tmp['log_C1_true']=-5
        tmp['log_C2_true']=-5
        tmp['r1_true']=100.*(.3)
        tmp['r2_true']=100.*(.3)
        tmp['log_alpha1_true'] =np.array(log_alpha)
        tmp['log_alpha2_true'] =np.array(log_alpha)
        tmp['log_gamma1_true'] =0.
        tmp['log_gamma2_true'] =0.
        T = T.append(tmp,ignore_index=True)
    except:
        print d

#Figure B, compare pso fitting
plt.figure(figsize=(3,5))
ax = []
cols = cm.jet(np.linspace(0,1,len(T['expt'].unique())))

fils = glob.glob('*/*/*.h5')
df = pd.DataFrame([])
x = {}
for tmp in T['expt'].unique():
    ky = tmp.split('_')[-1].split('x')[0]
    x[ky] = np.array([])
    x[ky+'_iter'] = np.array([])
    
    
for fil in fils:
    try:
        f = h5py.File(fil,'r')
        if 'pso' in f.keys():
            ky = fil.split('_')[1].split('x')[0]
            x[ky] = np.concatenate((x[ky],np.array((f['pso']['optimizer_logbook']['block0_values'])).T[0]))
            x[ky+'_iter'] = np.concatenate((x[ky+'_iter'],range(len(np.array((f['pso']['optimizer_logbook']['block0_values'])).T[0]))))
    except:
        print fil
        
        
num_expt = len(T['expt'].unique())
for e,tmp in enumerate(T['expt'].unique()):
    ax.append(plt.subplot2grid((num_expt,1),(e,0)))
    ky = tmp.split('_')[-1].split('x')[0]
    df = pd.DataFrame([x[ky],x[ky+'_iter']]).T

    y = df.groupby(1)[0].mean()
    y_std = df.groupby(1)[0].std()
    ax[e].plot(y.index.values,y.values,label=ky+'x'+ky,color=cols[e,:],linewidth=3)
    plt.sca(ax[e])
    plt.fill_between(y.index.values,y.values-y_std.values,y.values+y_std.values,alpha=0.5, color=cols[e,:],zorder=-int(ky))
    if e==0:
        ax[e].set_title('Log-likelihood of PSO-best')
    ax[e].set_ylabel(ky + 'x' + ky)
    ax[e].set_xlim([min(y.index.values),max(y.index.values)])
    ax[e].set_xticks([])
    ax[e].spines['top'].set_visible(False)
    ax[e].spines['right'].set_visible(False)
    ax[e].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    
ax[-1].set_xlabel('PSO Iteration')
ax[-1].set_xticks([min(y.index.values),max(y.index.values)])
plt.savefig('PSO-iter.pdf')
    

#Figure C, Compare the log-likelihood of different methods
fit_alg_order=['pso_mcmc','pso_multi_mcmc']
plt.figure(figsize=(5,2))
ax = [];
def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

for e,tmp in enumerate(T['expt'].unique()):
    ky = tmp.split('_')[-1].split('x')[0]
    ax.append(plt.subplot2grid((num_expt,1),(e,0)))
    sub_T = T[T['expt']==tmp]
    for e1,fa in enumerate(fit_alg_order):
        c = cols[e,:]
        c[-1]=(e1+1)/4. # change alpha
        if 'mcmc' in fa:
            indx = ((sub_T['fit_algorithm']==fa) & (sub_T['MCMC_converge']==1))
        else:
            indx = ((sub_T['fit_algorithm']==fa) & (sub_T['log_like_final']<100000000000.))

        data = sub_T['log_like_final'][indx].values
        if len(data)>5:
            parts = ax[e].violinplot(data,positions=[e1],showmeans=False, showmedians=False,showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(c)
                pc.set_edgecolor(c)
                pc.set_alpha(c[-1])
            quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=0)
            whiskers = adjacent_values(data, quartile1, quartile3)
            whiskersMin, whiskersMax = whiskers[0], whiskers[1]
            
            ax[e].scatter(e1, medians, marker='o', color='white', s=5, zorder=3)
            ax[e].vlines(e1, quartile1, quartile3, color='k', linestyle='-', lw=5)
            ax[e].vlines(e1, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1)
#            bp = ax[e].boxplot(data,notch=True,positions=[e1],showfliers=False)
#            for element in ['boxes', 'whiskers', 'means', 'medians', 'caps']:
#                plt.setp(bp[element], color=c)


    ax[e].set_ylabel(ky + 'x' + ky)
    ax[e].set_xlim([-.35,3.35])
    ax[e].set_xticks([])
    ax[e].spines['top'].set_visible(False)
    ax[e].spines['right'].set_visible(False)
    ax[e].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    
ax[-1].set_xticks(range(4))
ax[-1].set_xticklabels(fit_alg_order)
ax[0].set_title('R2 of fits')
plt.subplots_adjust(wspace=0,hspace=0.1)

    


#Figure D, Percent of models converged
plt.figure(figsize=(5,2))
ax = [];
for e,tmp in enumerate(T['expt'].unique()):
    ky = tmp.split('_')[-1].split('x')[0]
    ax.append(plt.subplot2grid((num_expt,1),(e,0)))
    sub_T = T[T['expt']==tmp]
    for e1,fa in enumerate(fit_alg_order[1:]):
        c = cols[e,:]
        c[-1]=(e1+1)/4. # change alpha
        data = sub_T[sub_T['fit_algorithm']==fa]
        if 'mcmc' in fa:
            val = float(sum((data['MCMC_converge']==1)&(data['R2']>.9)))/len(data)
            ax[e].bar(e1,val,color=c)
    if e == 0:
        ax[e].set_title('Percent Converged\n(Geweke Criterion)')
    
    ax[e].set_ylabel(ky + 'x' + ky)
    ax[e].set_ylim((.8,1.))
    ax[e].set_xlim([-.35,2.35])
    ax[e].set_xticks([])
    ax[e].spines['top'].set_visible(False)
    ax[e].spines['right'].set_visible(False)


labs = []
for fa in fit_alg_order:
    if 'mcmc' in fa:
        labs.append(fa)
ax[-1].set_xticks(range(3))
ax[-1].set_xticklabels(labs)






####Compare l2-norm for fitted parameter sets
fit_alg_order=['pso_mle','mcmc','pso_mcmc','pso_multi_mcmc']

nst_key = ['E0','E1','E2','E3','log_C1','log_C2','h1','h2','r1','r2','log_alpha1','log_alpha2']

indx = (T['MCMC_converge']==1) | ((T['fit_algorithm']=='pso_mle')&(T['log_like_final']<100000000000.))
sub_T = T[indx]
for k in nst_key:
    sub_T[k+'_diff'] = sub_T[k]-sub_T[k+'_true']
indx_key = list(sub_T.columns[sub_T.columns.str.contains('_diff')])
from sklearn import preprocessing
X_train = sub_T[indx_key]
X_scaled = preprocessing.scale(X_train)
dist = np.sqrt(np.sum(X_scaled**2,axis=1))
T['parm_distance']=np.nan
for e,ind in enumerate(sub_T.index):
    T.loc[ind,'parm_distance']=dist[e]
plt.figure(figsize=(3,5))
ax = [];
def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

for e,tmp in enumerate(T['expt'].unique()):
    ky = tmp.split('_')[-1].split('x')[0]
    ax.append(plt.subplot2grid((num_expt,1),(e,0)))
    sub_T = T[T['expt']==tmp]
    for e1,fa in enumerate(fit_alg_order):
        c = cols[e,:]
        c[-1]=(e1+1)/4. # change alpha
        if 'mcmc' in fa:
            if fa =='pso_multi_mcmc':
                indx = ((sub_T['fit_algorithm']==fa) & (sub_T['MCMC_converge']==1)&(sub_T['model_level']==5))
            else:
                indx = ((sub_T['fit_algorithm']==fa) & (sub_T['MCMC_converge']==1))
        else:
            indx = ((sub_T['fit_algorithm']==fa) & (sub_T['log_like_final']<100000000000.))

        data = sub_T['parm_distance'][indx].values
        if len(data)>5:
            parts = ax[e].violinplot(data,positions=[e1],showmeans=False, showmedians=False,showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(c)
                pc.set_edgecolor(c)
                pc.set_alpha(c[-1])
            quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=0)
            whiskers = adjacent_values(data, quartile1, quartile3)
            whiskersMin, whiskersMax = whiskers[0], whiskers[1]
            
            ax[e].scatter(e1, medians, marker='o', color='white', s=5, zorder=3)
            ax[e].vlines(e1, quartile1, quartile3, color='k', linestyle='-', lw=5)
            ax[e].vlines(e1, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1)
#            bp = ax[e].boxplot(data,notch=True,positions=[e1],showfliers=False)
#            for element in ['boxes', 'whiskers', 'means', 'medians', 'caps']:
#                plt.setp(bp[element], color=c)


    ax[e].set_ylabel(ky + 'x' + ky)
    ax[e].set_xlim([-.35,3.35])
    ax[e].set_xticks([])
    ax[e].spines['top'].set_visible(False)
    ax[e].spines['right'].set_visible(False)
    
ax[-1].set_xticks(range(4))
ax[-1].set_xticklabels(fit_alg_order)
ax[0].set_title('L2-norm parameter fits')
plt.subplots_adjust(wspace=0,hspace=0.2)
plt.savefig('L2-norm parameter fits.pdf')






######Plot DSD and where there are errors
cols2 = cm.jet(np.linspace(0,1,len(fit_alg_order)))
fit_alg_order=['pso_mcmc']

plt.figure(figsize=(5,2))
ax = [];
for e,tmp1 in enumerate(T['expt'].unique()):
    ky = tmp1.split('_')[-1].split('x')[0]
    plt.figure(figsize=(4,4))
    ax = plt.subplot(111)
    sub_T = T[T['expt']==tmp1]
    for e1,fa in enumerate(fit_alg_order):
        c = cols2[e1,:]
        if 'mcmc' in fa:
            indx = ((sub_T['fit_algorithm']==fa) & (sub_T['MCMC_converge']==1))
        else:
            indx = ((sub_T['fit_algorithm']==fa) & (sub_T['log_like_final']<100000000000.))

        data = sub_T[['beta_true','beta','log_alpha1','log_alpha1_true']][indx]
        ax.scatter(data['log_alpha1'],data['beta'],s=50,c=c,label=fa)
        for k in data.index:
            ax.plot(data[['log_alpha1','log_alpha1_true']].loc[k],data[['beta','beta_true']].loc[k],linewidth=.5,color='k',label='_nolegend_')
    ax.legend()
    tmp = np.unique(T[['log_alpha1_true','beta_true']],axis=0)
    ax.scatter(tmp[:,0],tmp[:,1],s=10,c='k',marker='+')
    ax.set_ylabel(r'$\beta$')
    ax.set_xlabel(r'$\alpha$')

    ax.set_xlim([-3,3])
    ax.set_ylim([-.7,.7])
    ax.set_title(ky + 'x' + ky)

plt.subplots_adjust(wspace=0,hspace=0.1)




##########################################################################
#Generate a heat map of alpha and beta uncertainty as a function of 
plt.figure(figsize=(4,4))
ax = plt.subplot(111)
for e,tmp1 in enumerate(T['expt'].unique()):
    ky = tmp1.split('_')[-1].split('x')[0]
    sub_T = T[(T['expt']==tmp1) & (T['fit_algorithm']=='pso_multi_mcmc') & (T['MCMC_converge']==1)]    
    y = sub_T.groupby('model_level')['model_level'].count()/len(sub_T)
    x1 = sum(y[y.index<5])
    x2 = 1-x1
    plt.bar(e,x1,color='r')
    plt.bar(e,x2,bottom=x1,color='b')
ax.set_xticks(range(len(T['expt'].unique())))   
ax.set_xticklabels(T['expt'].unique())   