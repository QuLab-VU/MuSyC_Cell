#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Code for estimating the % viability data from the HTS datasets
The cell count @ 0 and 72 hours are interpolated from the data.
Percent viability is estimated by normalizing
Cell Count @ (72hours) / Cell Count @ (0hours) * 2^(Basal proliferation rate * 72hours)
Basal proliferation rate is estimated from all the control wells.
Results are written to a csv file

Code commented and up to date 8/8/17
@author: Christian Meyer 
@email: christian.t.meyer@vanderbilt.edu
"""
##Import Packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from statsmodels import robust
import glob

fils = glob.glob('*.csv')
for f in fils:
    #Pick the time to measure the % viability
    measure_time = 72;
    expFile = f #Experiment file
    data = pd.read_csv(expFile)
    m_df = pd.DataFrame()
    for target_cell_line in np.unique(data['cell.line']):    
        #Read in data and subset to target cell line
        data = pd.read_csv(expFile)
        data = data[data['cell.line']==target_cell_line]
        
        #Identify the basal proliferation rate and standard deviation
        null_conc = np.mean(data[list(['drug1.conc','drug2.conc'])],axis=1) == 0
        df = data[null_conc]
        
        #For every control fit a line and estimate the slope and intercept.
        num_controls = len(np.unique(df['plate.name']))*len(np.unique(df['well']))
        slope = np.zeros((1,num_controls))
        inter = np.zeros((1,num_controls))
        sl_std = np.zeros((1,num_controls))
        inter_std = np.zeros((1,num_controls))
        cnt = 0
        #Plot the results
        plt.figure()
        #For each plate.
        for plate in np.unique(df['plate.name']):
            tmp = df[df['plate.name']==plate]
            tmp = tmp.reset_index()
            #For each well in each plate
            for well in np.unique(tmp['well']):
                tmp2 = tmp[tmp['well']==well];
                plt.plot(tmp2['time'],np.log2(tmp2['cell.count']))
                p,V = np.polyfit(tmp2['time'][tmp2['time']<40],np.log2(tmp2['cell.count'][tmp2['time']<40]),1,full = False, cov=True)
                slope[0,cnt] = p[0] 
                inter[0,cnt] = p[1] 
                sl_std[0,cnt] = V[0,0] 
                inter_std[0,cnt] = V[1,1]         
                cnt = cnt + 1
        
        #Allocate a data series for percent viability
        data['per_via'] = pd.Series(np.zeros((len(data),)), index=data.index)
        data['per_via_std'] = pd.Series(np.zeros((len(data),)), index=data.index)
    
        for plate in np.unique(data['plate.name']):
            tmp = data[data['plate.name']==plate]
            tmp = tmp.reset_index()
            for well in np.unique(tmp['well']):
                tmp2 = tmp[tmp['well']==well];
                plt.plot(tmp2['time'],np.log2(tmp2['cell.count']))
                s = interp1d(tmp2['time'], np.log2(tmp2['cell.count']), fill_value='extrapolate',kind='linear')
                zerotim = 2**(np.log2(tmp2[tmp2['time']==min(tmp2['time'])]['cell.count'])-np.median(slope)*min(tmp2['time']))
                zerotim_std = zerotim*robust.mad(slope,axis=1)/np.median(slope)
    #
                sv = 2**s(measure_time)
                c_tot = zerotim * 2**(np.median(slope)*measure_time)
                c_tot_std = np.sqrt((2**(np.median(slope)*measure_time)*zerotim_std)**2+
                                    (zerotim*np.log(2)*2**(np.median(slope)*measure_time)*robust.mad(slope,axis=1))**2)
                ind = data[np.logical_and(data['plate.name']==plate,data['well']==well)].index
                data.loc[ind,'per_via'] = float(sv/c_tot)
                data.loc[ind,'per_via_std'] = float(sv/c_tot*c_tot_std/c_tot)
        
        #Drop duplicate wells and plates.  This looses all the time information and just
        #keeps the first timepoint for every well.        
        df = data.drop_duplicates(subset=['plate.name','well'])
        
        #Remove some columns...
        #df.drop(list(['time','image.time','cell.count']), axis=1, inplace=True)
        #Write Results
        m_df =   pd.concat((m_df,df))  
        
    pd.DataFrame.to_csv(m_df,expFile + '_perVia_'+str(measure_time)+'hr.csv')
                
              