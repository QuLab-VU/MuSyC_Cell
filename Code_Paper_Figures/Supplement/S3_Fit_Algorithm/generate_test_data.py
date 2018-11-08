#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 13:43:32 2018

@author: xnmeyer
"""
import numpy as np
import pandas as pd
import os
import SynergyCalculator as SC
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
from shutil import copyfile 

#2D dose response surface function not assuming detail balance but assuming proliferation 
# is << than the rate of transition between the different states.
def Edrug2D_NDB(d,E0,E1,E2,E3,r1,r2,C1,C2,h1,h2,alpha1,alpha2):
    d1 = d[0]
    d2 = d[1]
    Ed = (E3*r1*(alpha1*d2)**h2*(alpha2*d1**2)**h1 + 
          E3*r2*(alpha2*d1)**h1*(alpha1*d2**2)**h2 + 
          C2**h2*E0*r1*(C1*alpha2*d1)**h1 + 
          C1**h1*E0*r2*(C2*alpha1*d2)**h2 + 
          C2**h2*E1*r1*(alpha2*d1**2)**h1 + 
          C1**h1*E2*r2*(alpha1*d2**2)**h2 + 
          E3*d2**h2*r1*(C1*alpha2*d1)**h1 + 
          E3*d1**h1*r2*(C2*alpha1*d2)**h2 + 
          C1**(2*h1)*C2**h2*E0*r1 + 
          C1**h1*C2**(2*h2)*E0*r2 + 
          C1**(2*h1)*E2*d2**h2*r1 + 
          C2**(2*h2)*E1*d1**h1*r2 + 
          E1*r2*(C2*d2)**h2*(alpha2*d1)**h1 + 
          E2*r1*(C1*d1)**h1*(alpha1*d2)**h2 +
          C2**h2*E1*r1*(C1*d1)**h1 +
          C1**h1*E2*r2*(C2*d2)**h2)/  \
         (r1*(alpha1*d2)**h2*(alpha2*d1**2)**h1 + 
          r2*(alpha2*d1)**h1*(alpha1*d2**2)**h2 + 
          C2**h2*r1*(C1*alpha2*d1)**h1 +
          C1**h1*r2*(C2*alpha1*d2)**h2 + 
          C2**h2*r1*(alpha2*d1**2)**h1 + 
          C1**h1*r2*(alpha1*d2**2)**h2 + 
          d2**h2*r1*(C1*alpha2*d1)**h1 + 
          d1**h1*r2*(C2*alpha1*d2)**h2 + 
          C1**(2*h1)*C2**h2*r1 + 
          C1**h1*C2**(2*h2)*r2 + 
          C1**(2*h1)*d2**h2*r1 + 
          C2**(2*h2)*d1**h1*r2 + 
          r1*(C1*d1)**h1*(alpha1*d2)**h2 + 
          r2*(C2*d2)**h2*(alpha2*d1)**h1 + 
          C2**h2*r1*(C1*d1)**h1 + 
          C1**h1*r2*(C2*d2)**h2)
         
    return Ed

exp_df = pd.read_csv('HTS018_rates.csv')
sd = exp_df['rate.95ci']/(2*1.96)
plt.figure()
ax = plt.subplot(111)
ax.hist(sd,bins=50)
noise = .001
data_mat = [5,7,10,15,25]
num_replicates = 2
for N in data_mat:
    df = pd.DataFrame()
    for e1,log_alpha in enumerate([-2,-1,0,1,2]):
        for e2,beta in enumerate([-.5,-.25,0,.25,.5]):
                for _ in range(num_replicates):
                    df_tmp = pd.DataFrame()
                    E0=.3
                    E1=0.
                    E2=0.
                    E3=-beta*(E0-min(E1,E2))+min(E1,E2) 
                    h1=1.
                    h2=1.
                    C1=10e-5
                    C2=10e-5
                    r1=100.*E0
                    r2=100.*E0
                    alpha1=10**log_alpha
                    alpha2=10**log_alpha
                    gamma1=1.
                    gamma2=1.
                    DD1,DD2 = np.meshgrid(np.concatenate(([0],np.logspace(-10,0,N))),np.concatenate(([0],np.logspace(-10,0,N))))
                    d1 = DD1.reshape((-1,))
                    d2 = DD2.reshape((-1,))
                    dip = Edrug2D_NDB((d1,d2),E0,E1,E2,E3,r1,r2,C1,C2,h1,h2,alpha1,alpha2)
                    dip = np.random.randn(len(dip))*noise+dip
                    dip_sd = np.ones(len(dip))*noise
                    df_tmp['drug1.conc']=d1
                    df_tmp['drug2.conc']=d2
                    df_tmp['effect.95ci'] =dip_sd
                    df_tmp['effect'] = dip
                    df_tmp['drug2'] = 'd2_alpha'+str(e1)+'_beta'+str(e2)
                    df_tmp['drug1'] = 'd1_alpha'+str(e1)+'_beta'+str(e2)
                    df_tmp['sample'] = 'NumData_' +str(N)+'_alpha'+str(e1)+'_beta'+str(e2)
                    df = df.append(df_tmp,ignore_index=True)
    
    os.mkdir('data_'+str(N)+'x'+str(N))
    for fit_alg in ['mcmc','pso_mcmc','multi_mcmc','pso_mle']: 
        os.mkdir('data_'+str(N)+'x'+str(N)+'/' + fit_alg)
        df.to_csv('data_'+str(N)+'x'+str(N)+'/' + fit_alg + '/' + 'testing_data_'+str(N)+'x'+str(N)+'.csv',index=False)
        

cnt = 1
for N in data_mat:
    for fit_alg in ['mcmc','pso_mcmc','multi_mcmc','pso_mle']:
        copyfile('SynergyCalculator.py','data_'+str(N)+'x'+str(N)+'/' + fit_alg + '/SynergyCalculator.py')
        copyfile('Syn_scratch.slurm','data_'+str(N)+'x'+str(N)+'/' + fit_alg + '/Syn_scratch.slurm')
        copyfile('SynCode_preCalcDIP_template.py','data_'+str(N)+'x'+str(N)+'/' + fit_alg + '/SynCode_preCalcDIP_template.py')

        os.chdir('data_'+str(N)+'x'+str(N)+'/' + fit_alg + '/')
        experiments = ['testing_data_'+str(N)+'x'+str(N)+'.csv']
        filename = 'SynCode_preCalcDIP_template'
        metric_name = 'DIP Rate'
        if fit_alg == 'multi_mcmc':
            min_nest  = 0
            fitting_alg = 'pso_mcmc'
        else:
            min_nest = 5
            fitting_alg = fit_alg
        max_nest = 5
        hill_orient = 1

        #Function to switch the index in SynCode.py
        def StringSwitch(string,file_string,val):
            idx1 = file_string.find(string)
            idx2 = file_string.find('#',idx1)
            file_string = file_string[0:idx1+len(string)] + str(val) + file_string[idx2:]
            return file_string    
    
        #Functions to generate array files by reading in the file name and saving an editted version with
        #a new index corresponding to a new drug combination.
        def Accre_Array_gen_darren(filename,target_cell_line,expt,drug1_name,drug2_name,metric_name,fitting_alg,cnt):
            fil = str(filename) + '.py'
            f = open(fil, 'r')
            file_string = f.read()
            string = 'fit_alg = '
            file_string = StringSwitch(string,file_string,'"%s"' % fitting_alg)
            string = 'expt = '
            file_string = StringSwitch(string,file_string,'"%s"' % expt)  
            string = 'metric_name = '
            file_string = StringSwitch(string,file_string,'"%s"' % metric_name)  
            string = 'min_nest = '
            file_string = StringSwitch(string,file_string,min_nest)
            string = 'max_nest = '
            file_string = StringSwitch(string,file_string,max_nest)
            string = 'hill_orient = '
            file_string = StringSwitch(string,file_string,hill_orient)  
            string = 'target_cell_line = '
            file_string = StringSwitch(string,file_string,'"%s"' % target_cell_line)
            string = 'drug1_name = '
            file_string = StringSwitch(string,file_string,'"%s"' % drug1_name)
            string = 'drug2_name = '
            file_string = StringSwitch(string,file_string,'"%s"' % drug2_name)  
            bol = True
            while bol:
                if os.path.isfile(filename + '_' + str(cnt) + '.py'):
                    cnt = cnt + 1
                else:
                    f = open(filename + '_' + str(cnt) + '.py', 'w')
                    f.write(file_string)
                    bol = False
            return cnt
        
        tmp_cnt = 0
        for expt in experiments:                
            #Read in the data into a pandas data frame
            data = pd.read_table(expt, delimiter=',')
            for e1,target_cell_line in enumerate(np.unique(data['sample'])):
                #Subset by target cell line
                sub_data = data[data['sample']==target_cell_line]
                drug_combinations = sub_data[list(['drug1','drug2'])].drop_duplicates()
                #Remove the double control condition...
                drug_combinations = drug_combinations[np.logical_and(drug_combinations['drug1']!='control', drug_combinations['drug2']!='control')]
                drug_combinations = drug_combinations.reset_index()
                #For all unique drug combinations:
                for e in drug_combinations.index:
                    drug1_name = drug_combinations['drug1'][e]
                    drug2_name = drug_combinations['drug2'][e]
                    cnt = Accre_Array_gen_darren(filename,target_cell_line,expt,drug1_name,drug2_name,metric_name,fitting_alg,cnt)          
                    if tmp_cnt==0:
                        init_cnt = cnt
                    tmp_cnt = 1

        fin_cnt = cnt
        cnt = cnt + 1
        f = open('Syn_scratch.slurm','r')
        file_string = f.read()
        idx1 = file_string.find('--array=')
        idx2 = file_string.find('\n',idx1)
        file_string = file_string[0:idx1+len('--array=')] + str(init_cnt) + '-' + str(fin_cnt) + file_string[idx2:]
        f.close()
        f = open('Syn_scratch.slurm','w')
        f.write(file_string)
        f.close()
        os.chdir('../..')
        
fil_str = 'cd mcmc/\nsbatch Syn_scratch.slurm\ncd ..\ncd multi_mcmc/\nsbatch Syn_scratch.slurm\ncd ..\ncd pso_mcmc/\nsbatch Syn_scratch.slurm\ncd ..\ncd pso_mle/\nsbatch Syn_scratch.slurm\ncd ..\n'
for N in data_mat:
    with open('submit_batch_files.txt','a') as f:
        f.write('cd data_'+str(N)+'x'+str(N)+'\n')
        f.write(fil_str)
        f.write('cd ..\n')
        

        
        
#############################################################################
### Now show how parameter uncertainty depends on dose response surface
#############################################################################

N  = 10
noise = .001
mx_dose = np.logspace(0,4,10)
for dose_surf in mx_dose:
    df = pd.DataFrame()
    for e,(log_alpha,beta) in enumerate(zip([-2,-1,1,2,0,0,0,0,0],[0,0,0,0,-.5,-.25,0,.25,.5])):
        df_tmp = pd.DataFrame()
        E0=.3
        E1=0.
        E2=0.
        E3=-beta*(E0-min(E1,E2))+min(E1,E2) 
        h1=1.
        h2=1.
        C1=10e-5
        C2=10e-5
        r1=100.*E0
        r2=100.*E0
        alpha1=10**log_alpha
        alpha2=10**log_alpha
        gamma1=1.
        gamma2=1.
        DD1,DD2 = np.meshgrid(np.concatenate(([0],np.logspace(-10,np.log10(dose_surf*C1),N))),np.concatenate(([0],np.logspace(-10,np.log10(dose_surf*C1),N))))
        d1 = DD1.reshape((-1,))
        d2 = DD2.reshape((-1,))
        dip = Edrug2D_NDB((d1,d2),E0,E1,E2,E3,r1,r2,C1,C2,h1,h2,alpha1,alpha2)
        dip = np.random.randn(len(dip))*noise+dip
        dip_sd = np.ones(len(dip))*noise
        df_tmp['drug1.conc']=d1
        df_tmp['drug2.conc']=d2
        df_tmp['effect.95ci'] =dip_sd
        df_tmp['effect'] = dip
        df_tmp['drug2'] = 'd2'+'_alpha'+str(log_alpha)+'_beta'+str(beta)
        df_tmp['drug1'] = 'd1'+'_alpha'+str(log_alpha)+'_beta'+str(beta)
        per_surf = np.log10(dose_surf*C1)
        df_tmp['sample'] = 'MaxDose_'+str(per_surf)+'_alpha'+str(log_alpha)+'_beta'+str(beta)
        df = df.append(df_tmp,ignore_index=True)
    
    os.mkdir('MaxDose_'+str(per_surf))
    for fit_alg in ['multi_mcmc','pso_mcmc']: 
        os.mkdir('MaxDose_'+str(per_surf)+'/' + fit_alg)
        df.to_csv('MaxDose_'+str(per_surf)+'/' + fit_alg + '/' + 'testing_maxdose_'+str(per_surf)+'.csv',index=False)
        

cnt = 1
for dose_surf in mx_dose:
    per_surf = np.log10(dose_surf*C1)
    for fit_alg in ['pso_mcmc']:
        copyfile('SynergyCalculator.py','MaxDose_'+str(per_surf)+'/' + fit_alg + '/SynergyCalculator.py')
        copyfile('Syn_scratch.slurm','MaxDose_'+str(per_surf)+'/' + fit_alg + '/Syn_scratch.slurm')
        copyfile('SynCode_preCalcDIP_template.py','MaxDose_'+str(per_surf)+'/' + fit_alg + '/SynCode_preCalcDIP_template.py')

        os.chdir('MaxDose_'+str(per_surf)+'/' + fit_alg + '/')
        experiments = ['testing_maxdose_'+str(per_surf)+'.csv']
        filename = 'SynCode_preCalcDIP_template'
        metric_name = 'DIP Rate'
        if fit_alg == 'multi_mcmc':
            min_nest  = 0
            fitting_alg = 'pso_mcmc'
        else:
            min_nest = 5
            fitting_alg = fit_alg
        max_nest = 5
        hill_orient = 1

        #Function to switch the index in SynCode.py
        def StringSwitch(string,file_string,val):
            idx1 = file_string.find(string)
            idx2 = file_string.find('#',idx1)
            file_string = file_string[0:idx1+len(string)] + str(val) + file_string[idx2:]
            return file_string    
    
        #Functions to generate array files by reading in the file name and saving an editted version with
        #a new index corresponding to a new drug combination.
        def Accre_Array_gen_darren(filename,target_cell_line,expt,drug1_name,drug2_name,metric_name,fitting_alg,cnt):
            fil = str(filename) + '.py'
            f = open(fil, 'r')
            file_string = f.read()
            string = 'fit_alg = '
            file_string = StringSwitch(string,file_string,'"%s"' % fitting_alg)
            string = 'expt = '
            file_string = StringSwitch(string,file_string,'"%s"' % expt)  
            string = 'metric_name = '
            file_string = StringSwitch(string,file_string,'"%s"' % metric_name)  
            string = 'min_nest = '
            file_string = StringSwitch(string,file_string,min_nest)
            string = 'max_nest = '
            file_string = StringSwitch(string,file_string,max_nest)
            string = 'hill_orient = '
            file_string = StringSwitch(string,file_string,hill_orient)  
            string = 'target_cell_line = '
            file_string = StringSwitch(string,file_string,'"%s"' % target_cell_line)
            string = 'drug1_name = '
            file_string = StringSwitch(string,file_string,'"%s"' % drug1_name)
            string = 'drug2_name = '
            file_string = StringSwitch(string,file_string,'"%s"' % drug2_name)  
            bol = True
            while bol:
                if os.path.isfile(filename + '_' + str(cnt) + '.py'):
                    cnt = cnt + 1
                else:
                    f = open(filename + '_' + str(cnt) + '.py', 'w')
                    f.write(file_string)
                    bol = False
            return cnt
        
        tmp_cnt = 0
        for expt in experiments:                
            #Read in the data into a pandas data frame
            data = pd.read_table(expt, delimiter=',')
            for e1,target_cell_line in enumerate(np.unique(data['sample'])):
                #Subset by target cell line
                sub_data = data[data['sample']==target_cell_line]
                drug_combinations = sub_data[list(['drug1','drug2'])].drop_duplicates()
                #Remove the double control condition...
                drug_combinations = drug_combinations[np.logical_and(drug_combinations['drug1']!='control', drug_combinations['drug2']!='control')]
                drug_combinations = drug_combinations.reset_index()
                #For all unique drug combinations:
                for e in drug_combinations.index:
                    drug1_name = drug_combinations['drug1'][e]
                    drug2_name = drug_combinations['drug2'][e]
                    cnt = Accre_Array_gen_darren(filename,target_cell_line,expt,drug1_name,drug2_name,metric_name,fitting_alg,cnt)          
                    if tmp_cnt==0:
                        init_cnt = cnt
                    tmp_cnt = 1

        fin_cnt = cnt
        cnt = cnt + 1
        f = open('Syn_scratch.slurm','r')
        file_string = f.read()
        idx1 = file_string.find('--array=')
        idx2 = file_string.find('\n',idx1)
        file_string = file_string[0:idx1+len('--array=')] + str(init_cnt) + '-' + str(fin_cnt) + file_string[idx2:]
        f.close()
        f = open('Syn_scratch.slurm','w')
        f.write(file_string)
        f.close()
        os.chdir('../..')
        
for dose_surf in mx_dose:
    per_surf = np.log10(dose_surf*C1)
    with open('submit_batch_files.txt','a') as f:
        f.write('cd MaxDose_'+str(per_surf)+'\n')
        f.write('cd pso_mcmc\nsbatch Syn_scratch.slurm\ncd ../..\n')
        




#
#import matplotlib.pylab as plt
#from mpl_toolkits.mplot3d import *
##from matplotlib.patches import FancyArrowPatch
#from matplotlib import rc
#rc('text', usetex=False)
#font = {'family' : 'arial',
#        'weight':'normal',
#        'size'   : 8}
#axes = {'linewidth': 2}
#rc('font', **font)
#rc('axes',**axes)
#
#fig = plt.figure()
#ax = plt.subplot(111,projection='3d')
#x = d1.copy()
#y = d2.copy()
#x[x==0]=min(x[x!=0])/10
#y[y==0]=min(y[y!=0])/10
#ax.scatter3D(np.log10(x),np.log10(y),dip)



