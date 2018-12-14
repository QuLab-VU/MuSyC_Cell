###
# HTS015 has very few points of osimertinib data, which makes it impossible to get accurate inference of the single-drug osimertinib dose-response
# Therefore, we will use the osi parameters from HTS018
# However, minor differences in the concentrations of osi between the HTS015 and HTS018 experiments mean that when we say conc=2e-7 in HTS015, it was more like 2e-6. 
# To compensate for this, and allow us to use osi parameters from HTS018 when calculating synergy in HTS015, we
# 1) In HTS015 data, for every concentration of osi alone, compute the median effect E_med
# 2) Using parameters from HTS018, and the inverse equation Edrug1D_inv, determine what the dose is that would give effect E_med
# 3) This dose becomes the dose we use when calculating Loewe and CI synergy of drugs in combination with osimertinib from HTS015
###

import pandas as pd
import numpy as np
import os


# 15: drug2 = osi
# 18: drug1 = osi

def Edrug1D(d,Em,C,h):
    return Em + (1-Em) / (1 + (d/C)**h)
    
def Edrug1D_inv(E, Em, C, h):
    if E > 1 or E < Em: return np.nan
    return C*np.power((1-Em)/(E-Em)-1., (1./h))


infile_directory = "../PerVia_Estimate"

df = pd.read_csv(os.path.join(infile_directory,"HTS015_dataset_PC9c1_sub.csv_perVia_72hr.csv"),index_col=0)
df['drug2.conc.original'] = df['drug2.conc']

params_18 = pd.read_csv("single_responses/18/parameters.csv",index_col=0)
params_15 = pd.read_csv("single_responses/15/parameters.csv",index_col=0)
osi_parameters = params_18.loc['osimertinib']
osi_params_df = pd.DataFrame(osi_parameters).transpose()
p15_index = list(params_15.index)
p15_index[-1] = "osimertinib_15"
params_15.index = p15_index
params_15 = pd.concat([params_15,osi_params_df])
params_15.to_csv("single_responses/15/parameters.csv")

    
drug = "osimertinib"
dff = df.loc[(df['drug2']==drug) & (df['drug1.conc']==0) & (df['drug2.conc']!=0),['drug2.conc','per_via']]

unique_concentrations = list(dff['drug2.conc'].unique())

Em = osi_parameters['emax']
C = osi_parameters['ec50']
h = osi_parameters['h']

transformed_concentrations = dict()

for d2c in unique_concentrations:
    E = dff.loc[dff['drug2.conc']==d2c,'per_via'].median()
    transformed_concentrations[d2c] = Edrug1D_inv(E, Em, C, h)
    df.loc[df['drug2.conc']==d2c,'drug2.conc'] = transformed_concentrations[d2c]
    
df.to_csv("results/HTS015_Per_viability_72hr_transformed.csv",index=False)
