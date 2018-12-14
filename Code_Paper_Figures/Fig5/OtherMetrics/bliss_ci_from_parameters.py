import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os

def Edrug1D(d,Em,C,h):
    return Em + (1-Em) / (1 + (d/C)**h)
    
def Edrug1D_inv(E, Em, C, h):
    if E > 1 or E < Em: return np.nan
    if np.isnan(C) or np.isnan(h): return np.nan
    return C*np.power((1-Em)/(E-Em)-1., (1./h))

def bliss_loewe(df,i,params,single_df):
    """
Bliss = 0 - Additive
Bliss > 0 - Synergistic
Bliss < 0 - Antagonistic 

Loewe = 1 - Additive
Loewe > 1 - Antagonistic
Loewe < 1 - Synergistic

CI = 1 - Additive
CI > 1 - Antagonistic
CI < 1 - Synergistic
"""
    drug1 = df.loc[i,'drug1']
    drug2 = df.loc[i,'drug2']
    d1c = df.loc[i,'drug1.conc']
    d2c = df.loc[i,'drug2.conc']
    
    plate = df.loc[i,'plate.id']
    
    d1_emax = params.loc[drug1,'emax']
    d1_ec50 = params.loc[drug1,'ec50']
    d1_h = params.loc[drug1,'h']
    
    d2_emax = params.loc[drug2,'emax']
    d2_ec50 = params.loc[drug2,'ec50']
    d2_h = params.loc[drug2,'h']
    
    
    d1_ec50_ci = params.loc[drug1,'ec50_CI']
    d1_h_ci = params.loc[drug1,'h_CI']


    d2_ec50_ci = params.loc[drug2,'ec50_CI']
    d2_h_ci = params.loc[drug2,'h_CI']
    
    combination = df.loc[i,'per_via']
    
    # Bliss (from curve)
#    e1_alone = Edrug1D(d1c,d1_emax,d1_ec50,d1_h)
#    e2_alone = Edrug1D(d2c,d2_emax,d2_ec50,d2_h)
#    bliss = e1_alone*e2_alone - combination
    
    # Bliss (directly from data)
    e1_alone = single_df.loc[(single_df['drug']==drug1) & (single_df['conc']==d1c) & (single_df['plate']==plate),'per_via'].mean()
    e2_alone = single_df.loc[(single_df['drug']==drug2) & (single_df['conc']==d2c) & (single_df['plate']==plate),'per_via'].mean()
    bliss = e1_alone*e2_alone - combination   # 0.9*0.9 - 0.5 ---> 0.81 - 0.5 = 0.31 ---->>> Synergistic
    
    # Loewe
    D1C = Edrug1D_inv(combination,d1_emax,d1_ec50,d1_h) # concentration of drug 1 alone which gets effect E
    D2C = Edrug1D_inv(combination,d2_emax,d2_ec50,d2_h) # concentration of drug 2 alone which gets effect E
    loewe = d1c/D1C + d2c/D2C
    
    # CI
    D1C = Edrug1D_inv(combination,0,d1_ec50_ci,d1_h_ci)
    D2C = Edrug1D_inv(combination,0,d2_ec50_ci,d2_h_ci)
    ci = d1c/D1C + d2c/D2C
    
    return bliss,loewe,ci
    
   


# HTS015
df = pd.read_csv("results/HTS015_Per_viability_72hr_transformed.csv")

single_drug_reference = df.loc[(df['drug1.conc']==0) ^ (df['drug2.conc']==0)] # '^' = exclusive or (XOR)
single_drug_dict = {'drug':[],'conc':[],'per_via':[], 'plate':[]}
for i in single_drug_reference.index:
    if single_drug_reference.loc[i,'drug1.conc']==0:
        single_drug_dict['drug'].append(single_drug_reference.loc[i,'drug2'])
        single_drug_dict['conc'].append(single_drug_reference.loc[i,'drug2.conc'])
    else:
        single_drug_dict['drug'].append(single_drug_reference.loc[i,'drug1'])
        single_drug_dict['conc'].append(single_drug_reference.loc[i,'drug1.conc'])
    single_drug_dict['per_via'].append(single_drug_reference.loc[i,'per_via'])
    single_drug_dict['plate'].append(single_drug_reference.loc[i,'plate.id'])

single_df = pd.DataFrame(single_drug_dict)

params = pd.read_csv("single_responses/15/parameters.csv",index_col=0)
df['bliss'] = np.nan
df['loewe_pervia'] = np.nan
df['CI'] = np.nan

for i in df.index:
    if df.loc[i,'drug1.conc']==0 or df.loc[i,'drug2.conc']==0: continue
    bliss,loewe,ci = bliss_loewe(df,i,params,single_df)
    df.loc[i,'bliss']=bliss
    df.loc[i,'loewe']=loewe
    df.loc[i,'CI']=ci

df['drug1'] = [i.lower() for i in df['drug1']]
df['drug2'] = [i.lower() for i in df['drug2']]
df.to_csv("results/HTS015_bliss_loewe_ci.csv",index=None)



#infile_directory = "/home/dwooten/School/qlab/Projects/dip-combination/dip_drug_synergy/data/HTS_Data"
infile_directory = "../PerVia_Estimate"
# HTS018
df = pd.read_csv(os.path.join(infile_directory,"HTS018_dataset_PC9c1_sub.csv_perVia_72hr.csv"),index_col=0)

single_drug_reference = df.loc[(df['drug1.conc']==0) ^ (df['drug2.conc']==0)] # '^' = exclusive or (XOR)
single_drug_dict = {'drug':[],'conc':[],'per_via':[],'plate':[]}
for i in single_drug_reference.index:
    if single_drug_reference.loc[i,'drug1.conc']==0:
        single_drug_dict['drug'].append(single_drug_reference.loc[i,'drug2'])
        single_drug_dict['conc'].append(single_drug_reference.loc[i,'drug2.conc'])
    else:
        single_drug_dict['drug'].append(single_drug_reference.loc[i,'drug1'])
        single_drug_dict['conc'].append(single_drug_reference.loc[i,'drug1.conc'])
    single_drug_dict['per_via'].append(single_drug_reference.loc[i,'per_via'])
    single_drug_dict['plate'].append(single_drug_reference.loc[i,'plate.id'])

single_df = pd.DataFrame(single_drug_dict)


params = pd.read_csv("single_responses/18/parameters.csv",index_col=0)
df['bliss'] = np.nan
df['loewe_pervia'] = np.nan
df['CI'] = np.nan

for i in df.index:
    if df.loc[i,'drug1.conc']==0 or df.loc[i,'drug2.conc']==0: continue
    bliss,loewe,ci = bliss_loewe(df,i,params,single_df)
    df.loc[i,'bliss']=bliss
    df.loc[i,'loewe']=loewe
    df.loc[i,'CI']=ci

df['drug1'] = [i.lower() for i in df['drug1']]
df['drug2'] = [i.lower() for i in df['drug2']]
df.to_csv("results/HTS018_bliss_loewe_ci.csv",index=None)




#infile_directory = "/home/dwooten/School/qlab/Projects/dip-combination/dip_drug_synergy/data/CellaVista_Data/per_via"

# Cellavista
df = pd.read_csv(os.path.join(infile_directory,"cellavista_pervia_combined.csv"))

single_drug_reference = df.loc[(df['drug1.conc']==0) ^ (df['drug2.conc']==0)] # '^' = exclusive or (XOR)
single_drug_dict = {'drug':[],'conc':[],'per_via':[], 'plate':[]}
for i in single_drug_reference.index:
    if single_drug_reference.loc[i,'drug1.conc']==0:
        single_drug_dict['drug'].append(single_drug_reference.loc[i,'drug2'])
        single_drug_dict['conc'].append(single_drug_reference.loc[i,'drug2.conc'])
    else:
        single_drug_dict['drug'].append(single_drug_reference.loc[i,'drug1'])
        single_drug_dict['conc'].append(single_drug_reference.loc[i,'drug1.conc'])
    single_drug_dict['per_via'].append(single_drug_reference.loc[i,'per_via'])
    single_drug_dict['plate'].append(single_drug_reference.loc[i,'plate.id'])

single_df = pd.DataFrame(single_drug_dict)

params = pd.read_csv("single_responses/cellavista/parameters.csv",index_col=0)
df['bliss'] = np.nan
df['loewe_pervia'] = np.nan
df['CI'] = np.nan

for i in df.index:
    if df.loc[i,'drug1.conc']==0 or df.loc[i,'drug2.conc']==0: continue
    bliss,loewe,ci = bliss_loewe(df,i,params,single_df)
    df.loc[i,'bliss']=bliss
    df.loc[i,'loewe']=loewe
    df.loc[i,'CI']=ci

df['drug1'] = [i.lower() for i in df['drug1']]
df['drug2'] = [i.lower() for i in df['drug2']]
df.to_csv("results/cellavista_bliss_loewe_ci.csv",index=None)
