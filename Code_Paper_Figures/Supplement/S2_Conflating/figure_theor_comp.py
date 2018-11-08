#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 11:13:10 2018

@author: xnmeyer
"""

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
from plotting_surfaces import plot_surface, plot_slices
from plotting_surfaces import get_params
from plotting_surfaces import dip
from pythontools.synergy import synergy_tools, synergy_plots
from scipy.interpolate import interp1d

fit_params = {}
fit_params['min_conc_d1'] = np.power(10.,-10)
fit_params['min_conc_d2'] = np.power(10.,-10)
fit_params['max_conc_d1'] = np.power(10.,0)
fit_params['max_conc_d2'] = np.power(10.,0)
fit_params['h1'] = .8
fit_params['h2'] = .8
fit_params['log_C1'] = -5
fit_params['log_C2'] = -5
fit_params['E0'] = 1
fit_params['E1'] =0.
fit_params['E2'] =0.
fit_params['E3'] = 0.0
fit_params['log_alpha1']=1.
fit_params['log_alpha2']=1.
fit_params['log_gamma1']=1.
fit_params['log_gamma2']=1.
fit_params['r1']=100.*fit_params['E0']
fit_params['r2']=100.*fit_params['E0']
fit_params['drug1_name']='d1'
fit_params['drug2_name']='d2'
fit_params['expt'] = 'cartoon'
fit_params = pd.DataFrame([fit_params])
N=100
zero_conc=1
zlim=(-0.05,1.05)
plot_surface(fit_params,incl_gamma=True,metric_name="Drug Effect",N=N,zero_conc=zero_conc,zlim=zlim,elev=17, azim=7)
plot_slices(fit_params,dd=1,zlim=zlim,incl_gamma=True,gN=N, N=10,fname=None, zero_conc=zero_conc,metric_name='%% Viable',title=None)  


E0, E1, E2, E3, h1, h2, r1, r1r, r2, r2r, C1, C2, alpha1, alpha2,gamma1,gamma2, d1min, d1max, d2min, d2max, drug1_name, drug2_name, expt = get_params(fit_params)

# Make plot showing Bliss and Loewe iso-contours in DSD space (calculate loewe and bliss at EC50)
alpha_range = np.logspace(-2,2,20)
beta_range = np.linspace(-.5,.5,20)

#Dose by dose
l_surf = np.zeros((len(beta_range),len(alpha_range)))
b_surf = np.zeros((len(beta_range),len(alpha_range)))
hsa_surf = np.zeros((len(beta_range),len(alpha_range)))
sch_surf = np.zeros((len(beta_range),len(alpha_range)))
#Whole surface
# 1D numpy arrays containing doses of drug1 and 2
d1 = np.concatenate((np.zeros(1),np.logspace(d1min,d1max,20))) # 10 doses
d2 = np.concatenate((np.zeros(1),np.logspace(d2min,d2max,20))) # 10 doses
DD1, DD2 = np.meshgrid(d1,d2)

brd_surf = np.zeros((len(beta_range),len(alpha_range)))
zip_surf = np.zeros((len(beta_range),len(alpha_range)))
zim_surf = np.zeros((len(beta_range),len(alpha_range)))

for i,alpha in enumerate(alpha_range):
    print i
    for j,beta in enumerate(beta_range):
        E3 = min(E1,E2) - beta*(E0-min(E1,E2))
        E = synergy_tools.hill_2D(C1, C2, E0, E1, E2, E3, h1, h2, alpha, alpha, r1, r1r, r2, r2r)
        l = -np.log10(synergy_tools.loewe(C1, C2, E, E0, E1, E2, h1, h2, C1, C2))
        l_surf[j,i] = l
        
        e_drug1_alone = synergy_tools.hill_1D(C1, E0, E1, h1, C1)
        e_drug2_alone = synergy_tools.hill_1D(C2, E0, E2, h2, C2)
        b = synergy_tools.bliss(e_drug1_alone, e_drug2_alone, E)
        b_surf[j,i] = b
        
        sch_surf[j,i] = synergy_tools.schindler(C1, C2, E, E0, E1, E2, h1, h2, C1, C2, logspace=False)
                
        E = synergy_tools.hill_2D(DD1, DD2, E0, E1, E2, E3, h1, h2, alpha, alpha, r1, r1r, r2, r2r)

        brd_surf[j,i] = synergy_tools.braid(DD1.reshape(-1,),DD2.reshape(-1,),E.reshape(-1,))
        hsa_surf[j,i] = synergy_tools.hsa(DD1.reshape(-1,),DD2.reshape(-1,),E.reshape(-1,))[0]
        zip_surf[j,i] = synergy_tools.zip(DD1,DD2,E)


import pickle
#with open('calculated_synergy_parms.p','w') as fil:
#    pickle.dump(brd_surf,fil)
#    pickle.dump(zip_surf,fil)
#    pickle.dump(hsa_surf,fil)
#    pickle.dump(sch_surf,fil)
#    pickle.dump(l_surf,fil)
#    pickle.dump(b_surf,fil)
#    
#to load back in
def loadall(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break
items = loadall('calculated_synergy_parms.p')
brd_surf = items.next()
zip_surf = items.next()
hsa_surf = items.next()
sch_surf = items.next()
l_surf = items.next()
b_surf = items.next()
    
ar = alpha_range;           br = beta_range
xlab = r"$\log(\alpha)$";   ylab = r"$\beta$"
CS,ax = synergy_plots.plot_surface_contour(ar, br, l_surf  , semilogx=True, xlabel=xlab, ylabel=ylab,crosshairs=True, levels=None, title="Loewe", power=2., figsize=(3.25,3.25), cmap ='plasma')
ax.set_yticks([-.5,0,.5])
ax.set_xticks([.01,1,100])
#A compicated way of identifying iso-surface locations
sav = 0
for e1,cl in enumerate(CS.collections):
    x = cl.get_paths()[0].vertices[:,0]
    y = cl.get_paths()[0].vertices[:,1]
    if x.min()<1 and x.max()>1 and y.min()<0 and y.max()>0:
        f1 = interp1d(x,y)
        f2 = interp1d(y,x)
        if sav<f2(0)/100+f1(1):
            sav = f2(0)/100+f1(1)
            sav_id = e1
x = CS.collections[sav_id].get_paths()[0].vertices[:,0]
y = CS.collections[sav_id].get_paths()[0].vertices[:,1]
f1 = interp1d(x,y);y1 = f1(1)
f2 = interp1d(y,x);x1 = f2(0)
lab1 = r'Loewe={0:.2f}'.format(CS.levels[sav_id]) + '\n' + r'$log(\alpha)$=0.00' + '\n' + r'$\beta$={0:.2f}'.format(float(y1))
lab2 = r'Loewe={0:.2f}'.format(CS.levels[sav_id]) + '\n' + r'$log(\alpha)$={0:.2f}'.format(float(np.log10(x1))) + '\n' + r'$\beta$=0.00'
ax.scatter(1,y1,s=50,c='r',marker='o',zorder=10,label=lab1)
ax.scatter(x1,0,s=50,c='r',marker='x',zorder=10,label=lab2)
plt.legend(loc='lower left')
plt.savefig('Plots/Loewe_contours.pdf',pad_inches=0.,format='pdf')

tmp = min([fit_params['E1'].values,fit_params['E2'].values])[0]
fit_params['E3'] = tmp
fit_params['log_alpha1']=np.log10(x1)
fit_params['log_alpha2']=np.log10(x1)
title = "Loewe:{0:.2f}".format(CS.levels[sav_id])+"\n"+r"$\beta$=0.0"+"\n"+r"$log(\alpha)$={0:.2f}".format(float(np.log10(x1)))
plot_surface(fit_params,metric_name="Drug Effect",N=N,zero_conc=zero_conc,zlim=(-.4,1.05),elev=41, azim=7,title=title,fname='Plots/loewe_alpha_surf')

tmp = min([fit_params['E1'].values,fit_params['E2'].values])[0]
fit_params['E3'] = tmp-y1*(fit_params['E0'].values[0]-tmp)
fit_params['log_alpha1']=0.
fit_params['log_alpha2']=0.
title = "Loewe:{0:.2f}".format(CS.levels[sav_id])+"\n"+r"$\beta$={0:.2f}".format(float(y1))+"\n"+r"$log(\alpha)$=0"
plot_surface(fit_params,metric_name="Drug Effect",N=N,zero_conc=zero_conc,zlim=(-.4,1.05),elev=41, azim=7,title=title,fname='Plots/loewe_beta_surf')

CS,ax = synergy_plots.plot_surface_contour(ar, br, sch_surf  , semilogx=True, xlabel=xlab, ylabel=ylab, crosshairs=True, levels=20, title="Schindler"  , power=1., figsize=(3.25,3.25), cmap ='plasma')
ax.set_yticks([-.5,0,.5])
ax.set_xticks([.01,1,100])
#A compicated way of identifying iso-surface locations
sav = 0
for e1,cl in enumerate(CS.collections):
    x = cl.get_paths()[0].vertices[:,0]
    y = cl.get_paths()[0].vertices[:,1]
    if x.min()<1 and x.max()>1 and y.min()<0 and y.max()>0:
        f1 = interp1d(x,y)
        f2 = interp1d(y,x)
        if sav<f2(0)/100+f1(1):
            sav = f2(0)/100+f1(1)
            sav_id = e1
x = CS.collections[sav_id].get_paths()[0].vertices[:,0]
y = CS.collections[sav_id].get_paths()[0].vertices[:,1]
f1 = interp1d(x,y);y1 = f1(1)
f2 = interp1d(y,x);x1 = f2(0)
lab1 = r'Schindler={0:.2f}'.format(CS.levels[sav_id]) + '\n' + r'$log(\alpha)$=0.00' + '\n' + r'$\beta$={0:.2f}'.format(float(y1))
lab2 = r'Schindler={0:.2f}'.format(CS.levels[sav_id]) + '\n' + r'$log(\alpha)$={0:.2f}'.format(float(np.log10(x1))) + '\n' + r'$\beta$=0.00'
ax.scatter(1,y1,s=50,c='r',marker='o',zorder=10,label=lab1)
ax.scatter(x1,0,s=50,c='r',marker='x',zorder=10,label=lab2)
plt.legend(loc='lower left')
plt.savefig('Plots/Schindler_contours.pdf',pad_inches=0.,format='pdf')

CS,ax = synergy_plots.plot_surface_contour(ar, br, hsa_surf  , semilogx=True, xlabel=xlab, ylabel=ylab, crosshairs=True, levels=20, title="HSA"  , power=1., figsize=(3.25,3.25), cmap ='plasma')
ax.set_yticks([-.5,0,.5])
ax.set_xticks([.01,1,100])
#A compicated way of identifying iso-surface locations
sav = 0
for e1,cl in enumerate(CS.collections):
    x = cl.get_paths()[0].vertices[:,0]
    y = cl.get_paths()[0].vertices[:,1]
    if x.min()<1 and x.max()>1 and y.min()<0 and y.max()>0:
        f1 = interp1d(x,y)
        f2 = interp1d(y,x)
        if sav<f2(0)/100+f1(1):
            sav = f2(0)/100+f1(1)
            sav_id = e1
x = CS.collections[sav_id].get_paths()[0].vertices[:,0]
y = CS.collections[sav_id].get_paths()[0].vertices[:,1]
f1 = interp1d(x,y);y1 = f1(1)
f2 = interp1d(y,x);x1 = f2(0)
lab1 = r'HSA={0:.2f}'.format(CS.levels[sav_id]) + '\n' + r'$log(\alpha)$=0.00' + '\n' + r'$\beta$={0:.2f}'.format(float(y1))
lab2 = r'HSA={0:.2f}'.format(CS.levels[sav_id]) + '\n' + r'$log(\alpha)$={0:.2f}'.format(float(np.log10(x1))) + '\n' + r'$\beta$=0.00'
ax.scatter(1,y1,s=50,c='r',marker='o',zorder=10,label=lab1)
ax.scatter(x1,0,s=50,c='r',marker='x',zorder=10,label=lab2)
plt.legend(loc='lower left')
plt.savefig('Plots/HSA_contours.pdf',pad_inches=0.,format='pdf')


CS,ax = synergy_plots.plot_surface_contour(ar, br, b_surf  , semilogx=True, xlabel=xlab, ylabel=ylab, crosshairs=True, levels=20, title="Bliss"  , power=1., figsize=(3.25,3.25), cmap ='plasma')
ax.set_yticks([-.5,0,.5])
ax.set_xticks([.01,1,100])
#A compicated way of identifying iso-surface locations
sav = 0
for e1,cl in enumerate(CS.collections):
    x = cl.get_paths()[0].vertices[:,0]
    y = cl.get_paths()[0].vertices[:,1]
    if x.min()<1 and x.max()>1 and y.min()<0 and y.max()>0:
        f1 = interp1d(x,y)
        f2 = interp1d(y,x)
        if sav<f2(0)/100+f1(1):
            sav = f2(0)/100+f1(1)
            sav_id = e1
x = CS.collections[sav_id].get_paths()[0].vertices[:,0]
y = CS.collections[sav_id].get_paths()[0].vertices[:,1]
f1 = interp1d(x,y);y1 = f1(1)
f2 = interp1d(y,x);x1 = f2(0)
lab1 = r'Bliss={0:.2f}'.format(CS.levels[sav_id]) + '\n' + r'$log(\alpha)$=0.00' + '\n' + r'$\beta$={0:.2f}'.format(float(y1))
lab2 = r'Bliss={0:.2f}'.format(CS.levels[sav_id]) + '\n' + r'$log(\alpha)$={0:.2f}'.format(float(np.log10(x1))) + '\n' + r'$\beta$=0.00'
ax.scatter(1,y1,s=50,c='r',marker='o',zorder=10,label=lab1)
ax.scatter(x1,0,s=50,c='r',marker='x',zorder=10,label=lab2)
plt.legend(loc='lower left')
plt.savefig('Plots/Bliss_contours.pdf',pad_inches=0.,format='pdf')


CS,ax = synergy_plots.plot_surface_contour(ar, br, zip_surf  , semilogx=True, xlabel=xlab, ylabel=ylab, crosshairs=True, levels=20, title="Zip"  , power=1., figsize=(3.25,3.25), cmap ='plasma')
ax.set_yticks([-.5,0,.5])
ax.set_xticks([.01,1,100])
#A compicated way of identifying iso-surface locations
sav_id = 8
x = CS.collections[sav_id].get_paths()[0].vertices[:,0]
y = CS.collections[sav_id].get_paths()[0].vertices[:,1]
f1 = interp1d(x,y);y1 = f1(.01)
f2 = interp1d(y,x);x1 = f2(0.01)
lab1 = r'Zip={0:.2f}'.format(CS.levels[sav_id]) + '\n' + r'$log(\alpha)$=-2.00' + '\n' + r'$\beta$={0:.2f}'.format(float(y1))
lab2 = r'Zip={0:.2f}'.format(CS.levels[sav_id]) + '\n' + r'$log(\alpha)$={0:.2f}'.format(float(np.log10(x1))) + '\n' + r'$\beta$=0.01'
ax.scatter(.01,y1,s=50,c='r',marker='o',zorder=10,label=lab1)
ax.scatter(x1,0.01,s=50,c='r',marker='x',zorder=10,label=lab2)
plt.legend(loc='lower left')
plt.savefig('Plots/Zip_contours.pdf',pad_inches=0.,format='pdf')




CS,ax = synergy_plots.plot_surface_contour(ar, br, brd_surf  , semilogx=True, xlabel=xlab, ylabel=ylab, crosshairs=True, levels=10, title="Braid"  , power=1., figsize=(3.25,3.25), cmap ='plasma')
ax.set_yticks([-.5,0,.5])
ax.set_xticks([.01,1,100])
#A compicated way of identifying iso-surface locations
sav = 0
xx,yy = np.meshgrid(ar,br)
df = pd.DataFrame([xx.reshape((-1,)),yy.reshape((-1,)),brd_surf.reshape((-1))]).T
a1,b1,c1 = df.loc[217,[0,1,2]]
a2,b2,c2 = df.loc[381,[0,1,2]]
lab1 = r'Braid={0:.2f}'.format(c1) + '\n' + r'$log(\alpha)$={0:.2f}'.format(np.log10(a1)) + '\n' + r'$\beta$={0:.2f}'.format(b1)
lab2 = r'Braid={0:.2f}'.format(c2) + '\n' + r'$log(\alpha)$={0:.2f}'.format(np.log10(a2)) + '\n' + r'$\beta$={0:.2f}'.format(b2)
ax.scatter(a1,b1,s=50,c='r',marker='o',zorder=10,label=lab1)
ax.scatter(a2,b2,s=50,c='r',marker='x',zorder=10,label=lab2)
plt.legend(loc='lower left')
plt.savefig('Plots/Braid_contours.pdf',pad_inches=0.,format='pdf')




######################################################
######################################################
#Show Zimmer and CI make bad fits for data which does not get to 0
def ll2(d,C,h):
    return 1/(1+(d/C)**h)
def ll4(d, E0, E1, h, C):
    return E1 + (E0-E1) / (1. + (d/C)**h)

alpha = 1
E0 = 1.
E3 = 0.
E1 = 0.5
E2 = 0.
h1 = .8
h2 = .8
r1r = 1e-5**h1
r2r = 1e-5**h2
r2 = 1
r1 = 1
d1min = -10
d2min = -10
d1max = 0
d2max = 0
d1 = np.concatenate((np.zeros(1),np.logspace(d1min,d1max,20))) # 10 doses
d2 = np.concatenate((np.zeros(1),np.logspace(d2min,d2max,20))) # 10 doses
DD1, DD2 = np.meshgrid(d1,d2)

E = synergy_tools.hill_2D(DD1, DD2, E0, E1, E2, E3, h1, h2, alpha, alpha, r1, r1r, r2, r2r)
C1_CI,C2_CI,h1_CI,h2_CI = synergy_tools.combination_index_fit(DD1.reshape((-1,)),DD2.reshape((-1,)),E.reshape((-1,)))
C1_zim, C2_zim, h1_zim, h2_zim = synergy_tools.zimmer(DD1,DD2,E)[0:4]
fig= plt.figure(figsize=(3.5,3))
ax = []
ax.append(plt.subplot2grid((2,2),(0,0)))
ax.append(plt.subplot2grid((2,2),(1,0)))
ax.append(plt.subplot2grid((2,2),(0,1)))
ax.append(plt.subplot2grid((2,2),(1,1)))
plt.sca(ax[0])
plt.plot(np.log10(np.unique(DD1)),ll2(np.unique(DD1),C1_CI,h1_CI),'r',label='CI Fit')
plt.plot(np.log10(np.unique(DD1)),ll4(np.unique(DD1),E0,E1,h1,C1),'b',label='True Curve')
plt.legend()
plt.ylabel('Drug Effect')
plt.xlabel('log([drug])')
plt.ylim([-.05,1.05])
plt.title('Combination Index')


plt.sca(ax[2])
plt.plot(np.log10(np.unique(DD1)),ll2(np.unique(DD1),C1_zim,h1_zim),'r',label='Dose Eqiv. Fit')
plt.plot(np.log10(np.unique(DD1)),ll4(np.unique(DD1),E0,E1,h1,C1),'b',label='True Curve')
plt.legend()
plt.xticks([])
plt.yticks([])
plt.ylim([-.05,1.05])
plt.title('Dose Equivalence Model')

from sklearn.metrics import r2_score
E1_vec = np.linspace(0.,.5,100)
R2_zim = np.zeros(len(E1_vec))
R2_ci = np.zeros(len(E1_vec))

for e,i in enumerate(E1_vec):
    E1 = i
    E = synergy_tools.hill_2D(DD1, DD2, E0, E1, E2, E3, h1, h2, alpha, alpha, r1, r1r, r2, r2r)
    C1_CI,C2_CI,h1_CI,h2_CI = synergy_tools.combination_index_fit(DD1.reshape((-1,)),DD2.reshape((-1,)),E.reshape((-1,)))
    C1_zim, C2_zim, h1_zim, h2_zim = synergy_tools.zimmer(DD1,DD2,E)[0:4]
    y_pred_CI =     ll2(np.unique(DD1),C1_CI,h1_CI)
    y_pred_zim =     ll2(np.unique(DD1),C1_zim,h1_zim)
    y_true = ll4(np.unique(DD1),E0,E1,h1,C1)
    R2_ci[e] = r2_score(y_true,y_pred_CI)
    R2_zim[e] = r2_score(y_true,y_pred_zim)
    
plt.sca(ax[1])
plt.plot(E1_vec,R2_ci)
plt.xlabel('True Emax')
plt.ylim([-.05,1.05])
plt.ylabel('R2 of CI fit')

plt.sca(ax[3])
plt.plot(E1_vec,R2_zim)
plt.xlabel('True Emax')
plt.ylabel('R2 of Dose Equiv. fit')
plt.ylim([-.05,1.05])


plt.tight_layout()
plt.savefig('Fitting_errors.pdf')



######################################################
#########Theoretical comparisons between different methods conflating potency and efficacy
#import pandas as pd
#import numpy as np
#import matplotlib.pylab as plt
##from matplotlib.patches import FancyArrowPatch
#from matplotlib import rc
#rc('text', usetex=False)
#font = {'family' : 'arial',
#        'weight':'normal',
#        'size'   : 8}
#axes = {'linewidth': 2}
#rc('font', **font)
#rc('axes',**axes)
#from plotting_surfaces import plot_surface
#from plotting_surfaces import get_params
#from plotting_surfaces import dip
#from pythontools.synergy import synergy_tools, synergy_plots
#from scipy.interpolate import interp1d
#
#fit_params = {}
#fit_params['min_conc_d1'] = np.power(10.,-10)
#fit_params['min_conc_d2'] = np.power(10.,-10)
#fit_params['max_conc_d1'] = np.power(10.,-5)
#fit_params['max_conc_d2'] = np.power(10.,-5)
#fit_params['h1'] = 1.
#fit_params['h2'] = 1.
#fit_params['log_C1'] = -7.5
#fit_params['log_C2'] = -7.5
#fit_params['E0'] = .05
#fit_params['E1'] = .005
#fit_params['E2'] = .005
#fit_params['E3'] = -.025
#fit_params['log_alpha1']=0.
#fit_params['log_alpha2']=0.
#fit_params['r1']=100000.*fit_params['E0']
#fit_params['r2']=100000.*fit_params['E0']
#fit_params['drug1_name']='d1'
#fit_params['drug2_name']='d2'
#fit_params['expt'] = 'cartoon'
#fit_params = pd.DataFrame([fit_params])
#N=50
#zero_conc=2
#zlim=(-0.055,0.055)
#plot_surface(fit_params,metric_name="Drug Effect",N=N,zero_conc=zero_conc,zlim=zlim,elev=41, azim=7)
#
#E0, E1, E2, E3, h1, h2, r1, r1r, r2, r2r, C1, C2, alpha1, alpha2, d1min, d1max, d2min, d2max, drug1_name, drug2_name, expt = get_params(fit_params)
#
## 1D numpy arrays containing doses of drug1 and 2
#d1 = np.concatenate((np.zeros(1),np.logspace(d1min,d1max,20))) # 10 doses
#d2 = np.concatenate((np.zeros(1),np.logspace(d2min,d2max,20))) # 10 doses
## Meshgrid of doses
#DD1, DD2 = np.meshgrid(d1,d2)
#
####
## PLOT VARIOUS HEATMAPS AND SURFACES
####
#
## MuSyC
#E = synergy_tools.hill_2D(DD1, DD2, E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)


