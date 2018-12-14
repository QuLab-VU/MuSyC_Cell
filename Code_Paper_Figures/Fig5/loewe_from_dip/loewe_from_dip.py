import numpy as np
import pandas as pd

fig3steps = """
For Loewe:
 Just use DIP
 Fails for synergistic efficacy: Vindesine
 Show a delta surface plot, with region showing where Loewe is totally undefined

For CI:
 Reduces to Loewe, BUT requires percent inhibition, from 0 to 1
 Show crappy fit
 18/carfilzomib
 18/gefitinib
 18/methotrexate ***
 18/ZM447439
 15/carmustine
 15/dronedarone
 15/GSK9948
 15/Vindesine
 15/Vinorelbine
 
 
For Bliss:
 Requires percent "affected", which is different than percent "inhibition"
 Calculate U, A1, A2, A12 = g(d1,d2)
 
 




Here we compare to traditional metrics. 
"""


def Edrug1D(d,Em,C,h):
    return Em + (1-Em) / (1 + (d/C)**h)
    
def Edrug1D_inv(E, Em, C, h):
    if E > 1 or E < Em: return np.nan
    if np.isnan(C) or np.isnan(h): return np.nan
    return C*np.power((1-Em)/(E-Em)-1., (1./h))

def f_loewe(df,i,h1,h2,E1,E2,C1,C2):
    """
Loewe = 1 - Additive
Loewe > 1 - Antagonistic
Loewe < 1 - Synergistic
"""
    d1c = df.loc[i,'drug1.conc']
    d2c = df.loc[i,'drug2.conc']
    
    combination = df.loc[i,'rate']
    if combination < min(E1,E2): return np.nan

    # Loewe
    D1C = Edrug1D_inv(combination,E1,C1,h1) # concentration of drug 1 alone which gets effect E
    D2C = Edrug1D_inv(combination,E2,C2,h2) # concentration of drug 2 alone which gets effect E
    loewe = d1c/D1C + d2c/D2C
    
    return loewe
    
master_drug_results = pd.read_csv("../../../Data/MasterResults_PC9_filtered.csv")


############ HTS015 ###############
df = pd.read_csv("../../../Data/HTS015_017_Combined.csv")
df = df.loc[df['cell.line']=="PC9c1"]
drug_id = 1
idx = [df.loc[i,"drug%d"%drug_id].lower() in master_drug_results['drug%d_name'%drug_id].values for i in df.index]
df = df.loc[idx]

df['loewe'] = np.nan

for i in df.index:
    if df.loc[i,'drug1.conc']==0 or df.loc[i,'drug2.conc']==0: continue
    drug = df.loc[i,'drug%d'%drug_id].lower()
    if not drug in master_drug_results['drug%d_name'%drug_id].values: continue 
    C1 = np.power(10.,np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"log_C1"]))
    C2 = np.power(10.,np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"log_C2"]))
    e1 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"E1"])
    e2 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"E2"])
    h1 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"h1"])
    h2 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"h2"])

    df.loc[i,'loewe'] = f_loewe(df,i,h1,h2,e1,e2,C1,C2)

df['drug1'] = [i.lower() for i in df['drug1']]
df['drug2'] = [i.lower() for i in df['drug2']]
df.to_csv("HTS015_loewe.csv",index=None)



########## HTS018 #############
df = pd.read_csv("../../../Data/HTS018_rates.csv")
df = df.loc[df['cell.line']=="PC9c1"]
drug_id = 2
idx = [df.loc[i,"drug%d"%drug_id].lower() in master_drug_results['drug%d_name'%drug_id].values for i in df.index]
df = df.loc[idx]

df['loewe'] = np.nan

for i in df.index:
    if df.loc[i,'drug1.conc']==0 or df.loc[i,'drug2.conc']==0: continue
    drug = df.loc[i,'drug%d'%drug_id].lower()
    if not drug in master_drug_results['drug%d_name'%drug_id].values: continue 
    C1 = np.power(10.,np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"log_C1"]))
    C2 = np.power(10.,np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"log_C2"]))
    e1 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"E1"])
    e2 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"E2"])
    h1 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"h1"])
    h2 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"h2"])

    df.loc[i,'loewe'] = f_loewe(df,i,h1,h2,e1,e2,C1,C2)

df['drug1'] = [i.lower() for i in df['drug1']]
df['drug2'] = [i.lower() for i in df['drug2']]
df.to_csv("HTS018_loewe.csv",index=None)




############ Cellavista ################
df = pd.read_csv("../../../Data/cellavista_combined.csv")

drug_id = 1

df['loewe'] = np.nan

for i in df.index:
    if df.loc[i,'drug1.conc']==0 or df.loc[i,'drug2.conc']==0: continue
    drug = df.loc[i,'drug%d'%drug_id].lower()
    if not drug in master_drug_results['drug%d_name'%drug_id].values: continue 
    C1 = np.power(10.,np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"log_C1"]))
    C2 = np.power(10.,np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"log_C2"]))
    e1 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"E1"])
    e2 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"E2"])
    h1 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"h1"])
    h2 = np.float(master_drug_results.loc[master_drug_results['drug%d_name'%drug_id]==drug,"h2"])

    df.loc[i,'loewe'] = f_loewe(df,i,h1,h2,e1,e2,C1,C2)

df['drug1'] = [i.lower() for i in df['drug1']]
df['drug2'] = [i.lower() for i in df['drug2']]
df.to_csv("cellavista_loewe.csv",index=None)
