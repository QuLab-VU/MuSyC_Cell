#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:23:03 2018

@author: xnmeyer
"""

#Figure Generation for Parameter Uncertainty
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
from scipy import stats


#Data generated in geneterate_test_data.py.  
#After fitting move to top directory
dirs = glob.glob('MaxDose*/*/')
#Compile all the results
T = pd.DataFrame([])
alpha_vals = [-2,-1,1,2,0,0,0,0,0]
beta_vals = [0,0,0,0,-.5,-.25,0,.25,.5]
for d in dirs:
    try:
        tmp = pd.read_csv(d+'MasterResults.csv')
        if 'DIC_0' in list(tmp) and tmp['fit_algorithm'].unique()[0]=='pso_mcmc':
            tmp['fit_algorithm']='pso_multi_mcmc'
        log_alpha = [float(i.split('_')[2][5:]) for i in tmp['cell_line']]
        beta = [float(i.split('_')[3][4:]) for i in tmp['cell_line']]
        tmp['beta_true'] = np.array(beta)
        tmp['E0_true'] = .3
        tmp['E1_true']=0.
        tmp['E2_true']=0.
        tmp['E3_true']=-np.array(beta)*.3 
        tmp['h1_true']=1.
        tmp['h2_true']=1.
        tmp['log_C1_true']=-4
        tmp['log_C2_true']=-4
        tmp['r1_true']=100.*(.3)
        tmp['r2_true']=100.*(.3)
        tmp['log_alpha1_true'] =np.array(log_alpha)
        tmp['log_alpha2_true'] =np.array(log_alpha)
        tmp['log_gamma1_true'] =0.
        tmp['log_gamma2_true'] =0.
        T = T.append(tmp,ignore_index=True)
    except:
        print d
        
    

#Calculate the uncertainty in log_alpha
fils = glob.glob('*/*/*.h5')
for fil in fils:
    arr = fil.split('_')
    mx_dose = float(arr[1].split('/')[0])
    fit_method = fil.split('/')[1]
    alpha = float(arr[5][5:])
    beta = float(arr[6][4:])
    indx = (T['fit_algorithm']==fit_method)&(T['beta_true']==beta)&(T['log_alpha1_true']==alpha)&(T['max_conc_d1']==10**mx_dose)    
    if sum(indx)==1:
        indx = T.index[np.where(indx)[0][0]]
        try:
            f = h5py.File(fil,'r')
            data = np.array(f['mcmc/trace_tier_5/block0_values'])
            keys = list(f['mcmc/trace_tier_5/block0_items'])
            for e,k in enumerate(keys):
                if k=='log_alpha1' or k=='log_alpha2':
                    x  = data[:,e]
                    T.loc[indx,k+'_std'] = np.std(data[1000:,e])
                    T.loc[indx,k]=np.mean(data[1000:,e])
        except IOError:
            continue


fils = glob.glob('*.0/*/*.h5')
#Pick out two example files
fils = fils[0],fils[14]
import pymc3 as pm
for fil in fils:
    arr = fil.split('_')
    mx_dose = float(arr[1].split('/')[0])
    fit_method = fil.split('/')[1]
    alpha = float(arr[5][5:])
    beta = float(arr[6][4:])
    indx = (T['fit_algorithm']==fit_method)&(T['beta_true']==beta)&(T['log_alpha1_true']==alpha)&(T['max_conc_d1']==10**mx_dose)    
    if sum(indx)==1:
        indx = T.index[np.where(indx)[0][0]]
        f = h5py.File(fil,'r')
        data = np.array(f['mcmc/trace_tier_5/block0_values'])
        keys = list(f['mcmc/trace_tier_5/block0_items'])
        plt.figure(figsize=(3.5,2))
        ax = []
        keys = list(f['mcmc/trace_tier_5/block0_items'])
        plt_keys = ['E3','log_alpha1','log_alpha2']
        for e,k in enumerate(keys):
            if np.in1d(k,plt_keys)[0]:
                ax.append(plt.subplot(len(plt_keys),3,3*plt_keys.index(k)+1))
                x  = data[:,e]
                ax[-1].plot(range(len(x)),x)
                ax[-1].set_ylabel(k)
                plt.axhline(T.loc[indx,k+'_true'],plt.xlim()[0],plt.xlim()[1],color='r',linestyle='--',label='True Value')
                if k.startswith('log_'):
                    plt.ylim((-.3,.1))
                else:
                    plt.ylim((-.2,.0))
                if k!='log_alpha2':
                    plt.xticks([])
                if k=='E3':
                    plt.title('MCMC Sampling')
                    plt.legend(loc='lower right')
                if k=='log_alpha2':
                    plt.xlabel('Iteration')
                    
                ax.append(plt.subplot(len(plt_keys),3,3*(plt_keys.index(k))+2))
                ge = pm.geweke(x,intervals=20)
                plt.scatter(ge[:,0],ge[:,1])
                x_lim = plt.xlim()
                plt.fill_between(x_lim,[-2,-2],[2,2],alpha=0.25, color='r',label=r'$Z-score\pm2$')
                plt.plot(x_lim,[2,2],'k-')
                plt.plot(x_lim,[-2,-2],'k-')
                if k=='E3':
                    plt.legend(loc='upper right')
                    plt.title('Geweke-score, interval=20')
                if k!='log_alpha2':
                    plt.xticks([])
                plt.ylim((-3,3))
                
                ax.append(plt.subplot(len(plt_keys),3,3*(plt_keys.index(k)+1)))
                x  = data[:,e]
                kernel = stats.gaussian_kde(x)
                if k.startswith('log_'):
                    plt.xlim((-.3,.1))
                else:
                    plt.xlim((-.2,.0))
                x = np.linspace(plt.xlim()[0],plt.xlim()[1],100)
                y = kernel(x)
                ax[-1].plot(x,y)
                ax[-1].set_ylabel(k)
                plt.axvline(T.loc[indx,k+'_true'],plt.ylim()[0],plt.ylim()[1],color='r',linestyle='--',label='True Value')
                if k!='log_alpha2':
                    plt.xticks([])
                if k=='E3':
                    plt.title('Posterior Distribution')

    plt.savefig('mx_dose_{:.2f}_traceplot.pdf'.format(mx_dose))
    
    
    


####Compare l2-norm for fitted parameter sets
plt.figure(figsize=(2,5))
ax = []
for em,key in enumerate(['log_alpha1','log_alpha2','beta']):
    fit_alg_order= ['pso_mcmc']
    ax.append(plt.subplot2grid((3,1),(em,0)))

    def adjacent_values(vals, q1, q3):
        upper_adjacent_value = q3 + (q3 - q1) * 1.5
        upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])
    
        lower_adjacent_value = q1 - (q3 - q1) * 1.5
        lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
        return lower_adjacent_value, upper_adjacent_value
    
    for e,tmp in enumerate(T['expt'].unique()):
        sub_T = T[T['expt']==tmp]
        for e1,fa in enumerate(fit_alg_order):
            x = np.log10(10**(float(tmp.split('_')[-1].split('.c')[0]))/(1e-5))
            if 'mcmc' in fa:
                if fa =='pso_multi_mcmc':
                    indx = ((sub_T['fit_algorithm']==fa) & (sub_T['MCMC_converge']==1)&(sub_T['R2']>.98))
                else:
                    indx = ((sub_T['fit_algorithm']==fa) & (sub_T['MCMC_converge']==1)&(sub_T['R2']>.98))
            else:
                indx = ((sub_T['fit_algorithm']==fa) & (sub_T['log_like_final']<100000000000.))
    
            data = (sub_T[key + '_std'][indx].values)
            if len(data)>5:
                quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=0)
                whiskers = adjacent_values(data, quartile1, quartile3)
                whiskersMin, whiskersMax = whiskers[0], whiskers[1]
                
                ax[em].scatter(e, medians, marker='o', color='white', s=5, zorder=3)
                ax[em].vlines(e, quartile1, quartile3, color='k', linestyle='-', lw=10)
                ax[em].vlines(e, whiskersMin, whiskersMax, color='k', linestyle='-', lw=2)

    
    ax[em].set_ylabel(key+'_std')
    ax[em].set_xlim([-.5,1.5])
    if em==2:
        ax[em].set_xticks(range(e+1))
        ax[em].set_xticklabels(['Max([d])=EC50','Max([d])=10000 x EC50'],rotation=45)
    else:
        ax[em].set_xticks([])
    ax[em].spines['top'].set_visible(False)
    ax[em].spines['right'].set_visible(False)
    plt.tight_layout()
plt.subplots_adjust(wspace=0,hspace=0.1)
plt.savefig('Uncertainty_dose_dependence.pdf')






#plot surface

########Theoretical comparisons between different methods conflating potency and efficacy
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



fils = glob.glob('*.0/*/*.h5')
#Pick out two example files
fils = fils[0],fils[14]
import pymc3 as pm
for fil in fils:
    arr = fil.split('_')
    mx_dose = float(arr[1].split('/')[0])
    fit_method = fil.split('/')[1]
    alpha = float(arr[5][5:])
    beta = float(arr[6][4:])
    indx = (T['fit_algorithm']==fit_method)&(T['beta_true']==beta)&(T['log_alpha1_true']==alpha)&(T['max_conc_d1']==10**mx_dose)    
    if sum(indx)==1:
        indx = T.index[np.where(indx)[0][0]]

    plot_surface(pd.DataFrame([T.loc[indx]]),data_pts=True,metric_name="DIP Rate",zlim=(-0.35,0.35),elev=24, azim=10,alpha=.9,fname=None,path=fil.split('2018')[0])
    plt.xlim((-11,0))
    plt.ylim((-11,0))
    plt.savefig('mx_dose_{:.2f}_surface.pdf'.format(mx_dose))













plt.savefig('L2-norm parameter fits.pdf')

        










cols = cm.hsv(np.linspace(0,1,len(T['beta_true'].unique())))

plt.figure()
ax = []
fit_alg = ['pso_mcmc']
for e,fa in enumerate(fit_alg):
    sub_T = T[T['MCMC_converge']==1 & (T['fit_algorithm']==fa) & (T['log_alpha1_true']==0.) & (T['R2']>.9)]
    ax.append(plt.subplot(1,2,e+1))
    x = sub_T.groupby(['beta_true','max_conc_d1'])[['beta_std']].mean().index.levels[0].values
    for e1,xx in enumerate(x):
        ax[-1].plot(sub_T.groupby(['beta_true','max_conc_d1'])[['beta_std']].mean().loc[xx],c=cols[e1],label=r'$\beta=$'+str(xx))
    ax[-1].legend()
    ax[-1].set_xscale('log')
    ax[-1].set_title(fa)
    
plt.figure()
ax = plt.subplot(121)
xx = .25
ax.plot(sub_T.groupby(['beta_true','max_conc_d1'])[['beta_std']].mean().loc[xx],c=cols[e1],label=r'$\beta=$'+str(xx))
ax.set_xscale('log')
  
ax = plt.subplot(122)
xx = 0.
ax.plot(sub_T.groupby(['log_alpha1_true','max_conc_d1'])[['alpha1_std']].mean().loc[xx],c=cols[e1],label=r'$\beta=$'+str(xx))
ax.set_xscale('log')

  
def conv_lognorm_to_norm_std(mu,sig):
    return -2*np.log10(mu)+np.log10(sig+mu)

def conv_norm_to_lognorm_std(mu,sig):
    return np.exp(2*mu+sig)*(np.exp(sig)-1)
    
T['log_alpha1_std'] = conv_norm_to_lognorm_std(10**T['log_alpha1'].values,T['alpha1_std'].values)
    
plt.figure()
ax = []
for e,fa in enumerate(fit_alg):
    sub_T = T[T['MCMC_converge']==1 & (T['fit_algorithm']==fa) & (T['beta_true']==0.) & (T['R2']>.9)]
    ax.append(plt.subplot(1,2,e+1))
    x = sub_T.groupby(['log_alpha1_true','max_conc_d1'])[['alpha1_std']].mean().index.levels[0].values
    for e1,xx in enumerate(x):
        ax[-1].plot(sub_T.groupby(['log_alpha1_true','max_conc_d1'])[['log_alpha1_std']].mean().loc[xx]/10**xx,c=cols[e1],label=r'$\log(alpha)=$'+str(xx))
    ax[-1].legend()
    ax[-1].set_xscale('log')
    ax[-1].set_title(fa)
    
    
fa = 'pso_mcmc'
sub_T = T[T['MCMC_converge']==1 & (T['fit_algorithm']==fa) & (T['log_alpha1_true']==0.) & (T['R2']>.7)]
ax.append(plt.subplot(1,2,e+1))
x = sub_T.groupby(['beta_true','max_conc_d1'])[['beta_std']].mean().index.levels[0].values
for e1,xx in enumerate(x):
    ax[-1].plot(sub_T.groupby(['beta_true','max_conc_d1'])[['beta_std']].mean().loc[xx],c=cols[e1],label=r'$\beta=$'+str(xx))
ax[-1].legend()
ax[-1].set_xscale('log')
ax[-1].set_title(fa)
    








fils = glob.glob('*/*/*.h5')
for fil in fils:
    f = h5py.File(fil,'r')
    data = np.array(f['mcmc/trace_tier_5/block0_values'])
    keys = list(f['mcmc/trace_tier_5/block0_items'])
    plt.figure(figsize=(5,10))
    ax = []
    keys = ['log_alpha1','E3','log_alpha2']
    for e,k in enumerate(keys):
        ax.append(plt.subplot(len(keys),2,2*e+1))
        x  = data[:,e]
        ax[-1].plot(range(len(x)),x)
        ax[-1].set_title(k)
        ax.append(plt.subplot(len(keys),2,2*(e+1)))
        ge = pm.geweke(x,intervals=20)
        plt.scatter(ge[:,0],ge[:,1])
        x_lim = plt.xlim()
        plt.plot(x_lim,[2,2],'k-')
        plt.plot(x_lim,[-2,-2],'k-')
        plt.title('Z-score for ' +  k)
    
    plt.tight_layout()
        








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


