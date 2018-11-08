# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 06:49:18 2018

@author: xnmeyer
"""
import pandas as pd
#Format tables for latex
t_pc9 = pd.read_csv('../../Data/MasterResults_PC9_filtered.csv')

keys = ['log_alpha1','log_alpha1_std','log_alpha2','log_alpha2_std','beta_obs_norm','beta_obs_norm_std','E1_obs','E1_obs_std','E2_obs','E2_obs_std','E3_obs','E3_obs_std']
for k in keys:
    t_pc9[k] = t_pc9[k].astype('float').round(3).astype('string')   

text = ''         
for i in t_pc9.index:
    if t_pc9['drug1_name'].loc[i]=='osimertinib':
        idx=2
        idt=1
    else:
        idx=1
        idt=2        
    text = text + \
           t_pc9['drug'+str(idx)+'_name'].loc[i]+'&$'+ \
           t_pc9['beta_obs_norm'].loc[i]+      '\pm'+t_pc9['beta_obs_norm_std'].loc[i]     +'$&$'+ \
           t_pc9['log_alpha'+str(idx)].loc[i]+ '\pm'+t_pc9['log_alpha'+str(idx)].loc[i]    +'$&$'+ \
           t_pc9['log_alpha'+str(idt)].loc[i]+ '\pm'+t_pc9['log_alpha'+str(idt)].loc[i]    +'$&$'+ \
           t_pc9['E'+str(idx)+'_obs'].loc[i]+  '\pm'+t_pc9['E'+str(idx)+'_obs_std'].loc[i] +'$&$'+ \
           t_pc9['E3_obs'].loc[i]+             '\pm'+t_pc9['E3_obs_std'].loc[i]            +'$\\\\'           
f = open('PC9_SuppTable7.txt','w')
f.write(text)
f.close()





T = pd.read_csv('../../Data/MasterResults.csv')
T = pd.read_csv('../../Data/MasterResults.csv')
T = T[T['cell_line']!='PC9C1']
T = T[T['cell_line']!='BR1']
T = T[T['cell_line']!='PC9AZR']
T = T[T['MCMC_converge']==1]
T = T[T['model_level']>4]

keys = ['log_alpha1','log_alpha1_std','log_alpha2','log_alpha2_std','beta_obs_norm','beta_obs_norm_std','E1_obs','E1_obs_std','E2_obs','E2_obs_std','E3_obs','E3_obs_std']
for k in keys:
    T[k] = T[k].astype('float').round(3).astype('string')   

T[r'$log(\alpha_1)$']   = '$' + T['log_alpha1'].astype('float').round(3).astype('string')       + '\pm' + T['log_alpha1_std'].astype('float').round(3).astype('string')     + '$'
T[r'$log(\alpha_2)$']   = '$' + T['log_alpha2'].astype('float').round(3).astype('string')       + '\pm' + T['log_alpha2_std'].astype('float').round(3).astype('string')     + '$'
T[r'$\beta_{obs}$']     = '$' + T['beta_obs_norm'].astype('float').round(3).astype('string')    + '\pm' + T['beta_obs_norm_std'].astype('float').round(3).astype('string')  + '$'
T[r'$E_1$']             = '$' + T['E1_obs'].astype('float').round(3).astype('string')           + '\pm' + T['E1_obs_std'].astype('float').round(3).astype('string')         + '$'
T[r'$E_2$']             = '$' + T['E2_obs'].astype('float').round(3).astype('string')           + '\pm' + T['E2_obs_std'].astype('float').round(3).astype('string')         + '$'
T[r'$E_3$']             = '$' + T['E3_obs'].astype('float').round(3).astype('string')           + '\pm' + T['E3_obs_std'].astype('float').round(3).astype('string')         + '$'
T['Cell-line']          =T['cell_line']
T['drug1']              =T['drug1_name']
T['drug2']              =T['drug2_name']

t = T.groupby(['Cell-line','drug1','drug2'])
f = open('MelanomaSynergyParameters_SuppTable6.txt','w')
f.write(t[[r'$\beta_{obs}$',r'$log(\alpha_1)$',r'$log(\alpha_2)$',r'$E_1$',r'$E_2$',r'$E_3$']].max().to_latex(escape=False))
f.close()



T = pd.read_csv('../../Data/MasterResults_PC9_filtered.csv')
############################################################################        
#Read in table of associated mechanism for each drug 
############################################################################
mech = pd.read_csv('../../Data/Drug_MOA_Table_supp.csv',index_col=0)
mech_key = pd.read_csv('../../Data/Drug_Class_Key.csv',index_col=0)
#Which ones to label in the DSDs
drug_combinations = t_pc9[['drug1_name','drug2_name','max_conc_d1','min_conc_d1','max_conc_d2','min_conc_d2']]      
#All drugs are in combination with osimertinib....Find the name of the other drug
drug_name=[]
tested_conc_rng=[]
for e,i in enumerate(drug_combinations.index):
    if drug_combinations['drug1_name'][i]!='osimertinib':
        drug_name.append(drug_combinations['drug1_name'][i])
        mx_conc = drug_combinations['max_conc_d1'][i]*1e6
        mn_conc = drug_combinations['min_conc_d1'][i]*1e9
        string = '{0:.1f}'.format(mx_conc) + 'uM-' + '{0:.1f}'.format(mn_conc) + 'nM,0nM'
        tested_conc_rng.append(string)
    else:
        drug_name.append(drug_combinations['drug2_name'][i])  
        mx_conc = drug_combinations['max_conc_d2'][i]*1e6
        mn_conc = drug_combinations['min_conc_d2'][i]*1e9
        string = '{0:.1f}'.format(mx_conc) + 'uM-' + '{0:.1f}'.format(mn_conc) + 'nM,0nM'
        tested_conc_rng.append(string)

#Name of super classes
sup_class = ['Mitotic\nCheckpoint', 'Epigenetic\nRegulators','Receptors&\nChannels', 'Kinases']
mech = mech.loc[drug_name]
mech['Drug Name']=mech.index


mech['Class'] = ''
mech['Subclass'] = ''
for i in mech_key.index:
    mech['Class'].loc[mech['label']==i]= mech_key['superclass'].loc[i]   
    mech['Subclass'].loc[mech['label']==i]=mech_key['subclass'].loc[i]

mech['Tested Range'] = tested_conc_rng

t = mech.groupby(['Class','Subclass','Drug Name','Tested Range','Nominal Target'])
f = open('PC9_drugPanel_Targets_suppTable5.txt','w')
f.write(t['Nominal Target'].all().to_latex(escape=False))
f.close()




##############################################################
##Read in the data and filter out the NSCLC cells to leave only the melanoma cells
###############################################################
T = pd.read_csv('../../Data/MasterResults.csv')
T = T[T['cell_line']!='PC9C1']
T = T[T['cell_line']!='BR1']
T = T[T['cell_line']!='PC9AZR']
T = T[T['MCMC_converge']==1]
T = T[T['model_level']>4]
T = T[T['R2']>.5]


t = T[['drug1_name','max_conc_d1','min_conc_d1']].groupby('drug1_name')
t.mean()

t = T[['drug2_name','max_conc_d2','min_conc_d2']].groupby('drug2_name')
t.mean()