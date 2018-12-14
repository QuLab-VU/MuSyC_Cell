#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 12:57:12 2018

@author: xnmeyer
"""
import pandas as pd

experiments = ['MasterResults.csv','MasterResults_plx_dpi_melPanel.csv']
for exp in experiments:
    T = pd.read_csv(exp)

    for ind in T.index:
        if ~np.isnan(T.loc[ind,['E1','E2']].values.astype('float')).all():
            #Determine which was the index of the most efficacious single agent
            minE = ['E1','E2'][(T.loc[ind,'E1'],T.loc[ind,'E2']).index(min(T.loc[ind,'E1'],T.loc[ind,'E2']))]
            #Define formula for calculating the standard deviation in beta
            def beta_std(T,minE):
                return np.sqrt(
                               ((T.loc[ind,'E0']-T.loc[ind,'E3'])/(T.loc[ind,'E0']-T.loc[ind,minE])**2*T.loc[ind,minE + '_std'])**2 + 
                               (T.loc[ind,'E3_std']/(T.loc[ind,'E0']-T.loc[ind,minE]))**2 + 
                               ((T.loc[ind,'E3']-T.loc[ind,minE])/(T.loc[ind,'E0']-T.loc[ind,minE])**2*T.loc[ind,'E0_std'])**2
                               )                
            
            T.loc[ind,'beta']                = (T.loc[ind,minE]-T.loc[ind,'E3'])/(T.loc[ind,'E0']-T.loc[ind,minE])
            T.loc[ind,'beta_std']            = beta_std(T,minE)
            #Single beta's
            T.loc[ind,'beta1']               = (T.loc[ind,'E1']-T.loc[ind,'E3'])/(T.loc[ind,'E0']-T.loc[ind,'E1'])
            T.loc[ind,'beta2']               = (T.loc[ind,'E2']-T.loc[ind,'E3'])/(T.loc[ind,'E0']-T.loc[ind,'E2'])
            T.loc[ind,'beta1_std']           = beta_std(T,'E1')
            T.loc[ind,'beta2_std']           = beta_std(T,'E2')
            
        if ~np.isnan(T.loc[ind,['E1_obs','E2_obs']].values.astype('float')).all():

            #Observed betas
            minE = ['E1_obs','E2_obs'][(T.loc[ind,'E1_obs'],T.loc[ind,'E2_obs']).index(min(T.loc[ind,'E1_obs'],T.loc[ind,'E2_obs']))]
            #Define formula for calculating the standard deviation in beta
            def beta_std(T,minE):
                return np.sqrt(
                               ((T.loc[ind,'E0']-T.loc[ind,'E3_obs'])/(T.loc[ind,'E0']-T.loc[ind,minE])**2*T.loc[ind,minE + '_std'])**2 + 
                               (T.loc[ind,'E3_obs_std']/(T.loc[ind,'E0']-T.loc[ind,minE]))**2 + 
                               ((T.loc[ind,'E3_obs']-T.loc[ind,minE])/(T.loc[ind,'E0']-T.loc[ind,minE])**2*T.loc[ind,'E0_std'])**2
                               )              
            
            T.loc[ind,'beta_obs']            = (T.loc[ind,minE]-T.loc[ind,'E3_obs'])/(T.loc[ind,'E0']-T.loc[ind,minE])
            T.loc[ind,'beta_obs_std']        = beta_std(T,minE)
            
            T.loc[ind,'beta1_obs']           = (T.loc[ind,'E1_obs']-T.loc[ind,'E3_obs'])/(T.loc[ind,'E0']-T.loc[ind,'E1_obs'])
            T.loc[ind,'beta2_obs']           = (T.loc[ind,'E2_obs']-T.loc[ind,'E3_obs'])/(T.loc[ind,'E0']-T.loc[ind,'E1_obs'])
            T.loc[ind,'beta1_obs_std']       = beta_std(T,'E1_obs')
            T.loc[ind,'beta2_obs_std']       = beta_std(T,'E2_obs')  
            
    T.to_csv(exp.split('.')[0],index=False)