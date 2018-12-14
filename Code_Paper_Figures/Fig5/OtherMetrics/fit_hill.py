import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
import os


# 15: drug2 = osi
# 18: drug1 = osi

# Assume that when d=0, E=1 (i.e., d=0 is control, so percent viable relative to control should be 1)
def Edrug1D(d,Em,C,h):
    return Em + (1-Em) / (1 + (d/C)**h)
    
def Edrug1D_log(d, Em, C, h):
    C_p = np.power(10.,C)
    d_p = np.power(10.,d)
    return Em + (1-Em) / (1 + (d_p/C_p)**h)


def fit_hill_log(d_conc, per_via, fname, plot=True):
    d_conc_log = np.log10(d_conc)
    try:
        popt1, pcov = curve_fit(Edrug1D_log, d_conc_log, per_via, bounds=([0,-np.inf,0],[1,0,np.inf]))
    
        if plot:
            dmin = d_conc_log.min()
            dmax = d_conc_log.max()
            plt.scatter(d_conc_log,per_via)
            x = np.linspace(dmin,dmax,50)
            plt.plot(x, Edrug1D_log(x, popt1[0], popt1[1], popt1[2]))
            plt.xlim(dmin,dmax)
            plt.ylim(0,max(1.1,1.1*per_via.max()))
            plt.savefig(fname)
            plt.close()
        popt1 = list(popt1)
        popt1[1] = np.power(10.,popt1[1])
        return list(popt1)
    except RuntimeError:
        print "Failed for drug %s"%drug
        dmin = d_conc_log.min()
        dmax = d_conc_log.max()
        plt.scatter(d_conc_log,per_via)
        plt.xlim(dmin,dmax)
        plt.ylim(0,max(1.1,1.1*per_via.max()))
        plt.title("FAILED")
        plt.savefig(fname)
        plt.close()
        return [np.nan, np.nan, np.nan]

def fit_hill_(d_conc, per_via, fname, plot=True):
    popt1, pcov = curve_fit(Edrug1D, d_conc, per_via, bounds=(0,[1,np.inf,np.inf]))

    if plot:
        dmin = np.log10(d_conc.min())
        dmax = np.log10(d_conc.max())
        plt.scatter(np.log10(d_conc),per_via)
        plt.plot(np.linspace(dmin,dmax,20), Edrug1D(np.logspace(dmin,dmax,20), popt1[0], popt1[1], popt1[2]))
        plt.xlim(dmin,dmax)
        plt.ylim(0,max(1.1,1.1*per_via.max()))
        plt.savefig(fname)
        plt.close()
    
    return list(popt1)
    

# Use scipy.optimize.curve_fit to fit sigmoidal Hill equation
def fit_hill_new(d_conc, per_via, fname, plot=False):
    try:
        if drug=="Vinorelbine": popt1, pcov = curve_fit(Edrug1D_log, np.log10(d_conc), per_via, bounds=([0,-7.75,0],[0.4,-7.,np.inf]))
        else: popt1, pcov = curve_fit(Edrug1D_log, np.log10(d_conc), per_via, bounds=([0,-np.inf,0],[1,0,np.inf]))
    except RuntimeError:
        print "\n\n\n*********\nfailed for drug %s\n*********\n\n\n\n"%drug
        if plot:
            dmin = np.log10(d_conc.min())
            dmax = np.log10(d_conc.max())
            plt.scatter(np.log10(d_conc),per_via)
            plt.xlim(dmin,dmax)
            plt.ylim(0,max(1.1,1.1*per_via.max()))
            plt.savefig(fname)
            plt.close()
        return [np.nan, np.nan, np.nan]
    if plot:
        dmin = np.log10(d_conc.min())
        dmax = np.log10(d_conc.max())
        plt.scatter(np.log10(d_conc),per_via)
        x = np.linspace(dmin,dmax,20)
        plt.plot(x, Edrug1D_log(x, popt1[0], popt1[1], popt1[2]))
        plt.xlim(dmin,dmax)
        plt.ylim(0,max(1.1,1.1*per_via.max()))
        plt.savefig(fname)
        plt.close()
    
    return list(popt1)

# CI assumes that log(fA / fU) is a linear function of log(dose), and only works for values with 0 <= fU,fA <= 1
def fit_hill_ci(d_conc, per_via, fname, plot=True):
    fU = per_via
    fA = 1-per_via
    C2 = linregress(np.log(d_conc),np.log(fA/fU))
    popt1 = [0,0,0]
    popt1[2] = C2.slope
    popt1[1] = np.exp(C2.intercept/(-popt1[2])) # IC50 for drug 1
    
    if plot:
        dmin = np.log10(d_conc.min())
        dmax = np.log10(d_conc.max())
        plt.scatter(np.log10(d_conc),per_via)
        plt.plot(np.linspace(dmin,dmax,20), Edrug1D(np.logspace(dmin,dmax,20), popt1[0], popt1[1], popt1[2]))
        plt.xlim(dmin,dmax)
        plt.ylim(0,max(1.1,1.1*per_via.max()))
        plt.savefig(fname)
        plt.close()
    
    return popt1


infile_directory = "../PerVia_Estimate"

# Fit hill curves to single-dose-responses 

# HTS018
df = pd.read_csv(os.path.join(infile_directory,"HTS018_dataset_PC9c1_sub.csv_perVia_72hr.csv"))
df = df.loc[df['cell.line']=="PC9c1"]
single_df = pd.DataFrame(index=['emax','ec50','h'])
df_ci = df.loc[df['per_via']<=1.]
single_df_ci = pd.DataFrame(index=['emax','ec50','h'])

for drug in df['drug2'].unique():
    if drug == "control": continue

    # Fit standard Hill
    dff = df.loc[(df['drug2']==drug) & (df['drug1.conc']==0),['drug2.conc','per_via']].sort_values(by="drug2.conc",ascending=True)
    single_df[drug]=fit_hill_log(dff['drug2.conc'], dff['per_via'], "single_responses/18/%s.png"%drug)
    
    # Fit Hill for CI
    dff = df_ci.loc[(df_ci['drug2']==drug) & (df_ci['drug1.conc']==0),['drug2.conc','per_via']].sort_values(by="drug2.conc",ascending=True)
    single_df_ci[drug]=fit_hill_ci(dff['drug2.conc'], dff['per_via'], "single_responses/CI/18/%s.png"%drug)


drug = "osimertinib"
# Fit standard Hill
dff = df.loc[(df['drug1']==drug) & (df['drug2']=="control"),['drug1.conc','per_via']]
single_df[drug]=fit_hill_log(dff['drug1.conc'], dff['per_via'], "single_responses/18/%s.png"%drug)

# Fit Hill for CI
dff = df_ci.loc[(df_ci['drug1']==drug) & (df_ci['drug2']=="control"),['drug1.conc','per_via']]
single_df_ci[drug]=fit_hill_ci(dff['drug1.conc'], dff['per_via'], "single_responses/CI/18/%s.png"%drug)


single_df = single_df.transpose()
#single_df['ec50_p'] = np.power(10.,single_df['ec50'])

single_df_ci = single_df_ci.transpose()
single_df['ec50_CI'] = single_df_ci['ec50']
single_df['h_CI'] = single_df_ci['h']
single_df.to_csv("single_responses/18/parameters.csv")
single_df_ci.to_csv("single_responses/CI/18/parameters.csv")



# HTS015
df = pd.read_csv(os.path.join(infile_directory,"HTS015_dataset_PC9c1_sub.csv_perVia_72hr.csv"),index_col=0)
df = df.loc[df['cell.line']=="PC9c1"]
single_df = pd.DataFrame(index=['emax','ec50','h'])
df_ci = df.loc[df['per_via']<=1.]
single_df_ci = pd.DataFrame(index=['emax','ec50','h'])
for drug in df['drug1'].unique():
    if drug=="control": continue

    # Fit standard Hill
    dff = df.loc[(df['drug1']==drug) & (df['drug1.conc']!=0) & (df['drug2.conc']==0),['drug1.conc','per_via']].sort_values(by="drug1.conc",ascending=True)
    single_df[drug]=fit_hill_log(dff['drug1.conc'], dff['per_via'], "single_responses/15/%s.png"%drug)
    
    # Fit CI Hill
    dff = df_ci.loc[(df_ci['drug1']==drug) & (df_ci['drug1.conc']!=0) & (df_ci['drug2.conc']==0),['drug1.conc','per_via']].sort_values(by="drug1.conc",ascending=True)
    try:
        single_df_ci[drug]=fit_hill_ci(dff['drug1.conc'], dff['per_via'], "single_responses/CI/15/%s.png"%drug)
    except ValueError:
        single_df_ci[drug] = [0,np.nan,np.nan]

drug = "osimertinib"

#Fit standard Hill
dff = df.loc[(df['drug2']==drug) & (df['drug1.conc']==0) & (df['drug2.conc']!=0),['drug2.conc','per_via']]
single_df[drug]=fit_hill_log(dff['drug2.conc'], dff['per_via'], "single_responses/15/%s.png"%drug)

#Fit CI Hill
dff = df_ci.loc[(df_ci['drug2']==drug) & (df_ci['drug1.conc']==0) & (df_ci['drug2.conc']!=0),['drug2.conc','per_via']]
try:
    single_df_ci[drug]=fit_hill_ci(dff['drug2.conc'], dff['per_via'], "single_responses/CI/15/%s.png"%drug)
except ValueError:
    single_df_ci[drug]=[0,np.nan,np.nan]

single_df = single_df.transpose()
#single_df['ec50'] = np.power(10.,single_df['ec50'])
single_df_ci = single_df_ci.transpose()
single_df['ec50_CI'] = single_df_ci['ec50']
single_df['h_CI'] = single_df_ci['h']
single_df.to_csv("single_responses/15/parameters.csv")
single_df_ci.to_csv("single_responses/CI/15/parameters.csv")
    
   
   

   
# Cellavista
df = pd.read_csv(os.path.join(infile_directory,"cellavista_pervia_combined.csv"),index_col=0)
single_df = pd.DataFrame(index=['emax','ec50','h'])
df_ci = df.loc[df['per_via']<=1.]
single_df_ci = pd.DataFrame(index=['emax','ec50','h'])
for drug in df['drug1'].unique():
    if drug=="control": continue
    
    # Fit standard Hill
    dff = df.loc[(df['drug1']==drug) & (df['drug1.conc']!=0) & (df['drug2.conc']==0),['drug1.conc','per_via']].sort_values(by="drug1.conc",ascending=True)
    single_df[drug]=fit_hill_log(dff['drug1.conc'], dff['per_via'], "single_responses/cellavista/%s.png"%drug)
    
    # Fit CI Hill
    dff = df_ci.loc[(df_ci['drug1']==drug) & (df_ci['drug1.conc']!=0) & (df_ci['drug2.conc']==0),['drug1.conc','per_via']].sort_values(by="drug1.conc",ascending=True)
    single_df_ci[drug]=fit_hill_ci(dff['drug1.conc'], dff['per_via'], "single_responses/CI/cellavista/%s.png"%drug)
    
drug = "osimertinib"

# Fit standard Hill
dff = df.loc[(df['drug2']==drug) & (df['drug1.conc']==0) & (df['drug2.conc']!=0),['drug2.conc','per_via']]
single_df[drug]=fit_hill_log(dff['drug2.conc'], dff['per_via'], "single_responses/cellavista/%s.png"%drug)

# Fit CI Hill
dff = df_ci.loc[(df_ci['drug2']==drug) & (df_ci['drug1.conc']==0) & (df_ci['drug2.conc']!=0),['drug2.conc','per_via']]
single_df_ci[drug]=fit_hill_ci(dff['drug2.conc'], dff['per_via'], "single_responses/CI/cellavista/%s.png"%drug)

single_df = single_df.transpose()
#single_df['ec50'] = np.power(10.,single_df['ec50'])
single_df_ci = single_df_ci.transpose()
single_df['ec50_CI'] = single_df_ci['ec50']
single_df['h_CI'] = single_df_ci['h']
single_df.to_csv("single_responses/cellavista/parameters.csv")
single_df_ci.to_csv("single_responses/CI/cellavista/parameters.csv")
