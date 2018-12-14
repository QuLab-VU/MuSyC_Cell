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


def fit_hill_log(d_conc, per_via, fname, plot=False):
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
def fit_hill_ci(d_conc, per_via, fname, plot=False):
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


def manuscript_plot(single_df, df, df_ci, drug, ax1, ax2, left):
#    fig = plt.figure()
    fU = df_ci['per_via']
    fA = 1-df_ci['per_via']
    x = np.log(df_ci['drug2.conc'])
    y = np.log(fA/fU)
    
    C2 = linregress(x,y)
    h_ci = C2.slope
    c_ci = np.exp(C2.intercept/(-h_ci))
    
    print h_ci, c_ci
    
    ax = ax1
    ax.scatter(x,y,c='k',s=20)
    ax.set_xlim(x.min(),x.max())
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.plot(xlim, C2.intercept + np.asarray(xlim)*C2.slope)
    if left: ax.set_ylabel(r"$\ln\left(\frac{f_A}{f_U}\right)$", fontsize=8)
#    ax.set_xlabel("Concentration\n" + r"($\ln[M]$)", fontsize=8)
    ax.set_xlim(xlim)
    ax.set_title(drug, fontsize=8)
    ax.tick_params(labelsize=8)
    
    
    
    ax = ax2
    
    d_conc_log = np.log10(df['drug2.conc'])

    dmin = d_conc_log.min()
    dmax = d_conc_log.max()
    dmax += 2.
    
    
    
    x = np.linspace(dmin,dmax,50)
    ax.plot(x, Edrug1D_log(x, 0., np.log10(float(single_df['ec50_ci'])), float(single_df['h_ci'])), label="CI")
    ax.plot(x, Edrug1D_log(x, float(single_df['emax']), np.log10(float(single_df['ec50'])), float(single_df['h'])), label="curve_fit")
#    ax.legend(fontsize=8)
    ax.scatter(d_conc_log,df['per_via'], c='k', s=20)
    ax.set_xlim(dmin,dmax)
    ax.set_ylim(0,max(1.1,1.1*df['per_via'].max()))

    if left: ax.set_ylabel("Percent Viability", fontsize=8)
#    ax.set_xlabel("Concentration\n" + r"($\log_{10}[M]$)", fontsize=8)
    ax.tick_params(labelsize=8)
    
    
#    plt.tight_layout()
#    plt.show()


infile_directory = "../PerVia_Estimate"

# Fit hill curves to single-dose-responses 

# HTS018
df = pd.read_csv(os.path.join(infile_directory,"HTS018_dataset_PC9c1_sub.csv_perVia_72hr.csv"))
df = df.loc[df['cell.line']=="PC9c1"]

single_df = pd.DataFrame(index=['emax','ec50','h'])
df_ci = df.loc[df['per_via']<=1.]
single_df_ci = pd.DataFrame(index=['emax','ec50','h'])



fig = plt.figure(figsize=(3,3))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)


drug = "panobinostat"
# Fit Hill for CI
dff_ci = df_ci.loc[(df_ci['drug2']==drug) & (df_ci['drug1.conc']==0),['drug2.conc','per_via']].sort_values(by="drug2.conc",ascending=True)
single_df_ci[drug]=fit_hill_ci(dff_ci['drug2.conc'], dff_ci['per_via'], "manuscript/%s_ci.png"%drug)

# Fit standard Hill
dff = df.loc[(df['drug2']==drug) & (df['drug1.conc']==0),['drug2.conc','per_via']].sort_values(by="drug2.conc",ascending=True)
single_df[drug]=fit_hill_log(dff['drug2.conc'], dff['per_via'], "manuscript/%s.png"%drug)


single_df = single_df.transpose()
single_df_ci = single_df_ci.transpose()

single_df['ec50_ci'] = single_df_ci['ec50']
single_df['h_ci'] = single_df_ci['h']
manuscript_plot(single_df, dff, dff_ci, drug, ax1, ax3, True)






drug = "methotrexate"
single_df = pd.DataFrame(index=['emax','ec50','h'])
df_ci = df.loc[df['per_via']<=1.]
single_df_ci = pd.DataFrame(index=['emax','ec50','h'])
# Fit Hill for CI
dff_ci = df_ci.loc[(df_ci['drug2']==drug) & (df_ci['drug1.conc']==0),['drug2.conc','per_via']].sort_values(by="drug2.conc",ascending=True)
single_df_ci[drug]=fit_hill_ci(dff_ci['drug2.conc'], dff_ci['per_via'], "manuscript/%s_ci.png"%drug)

# Fit standard Hill
dff = df.loc[(df['drug2']==drug) & (df['drug1.conc']==0),['drug2.conc','per_via']].sort_values(by="drug2.conc",ascending=True)
single_df[drug]=fit_hill_log(dff['drug2.conc'], dff['per_via'], "manuscript/%s.png"%drug)

single_df = single_df.transpose()
single_df_ci = single_df_ci.transpose()

single_df['ec50_ci'] = single_df_ci['ec50']
single_df['h_ci'] = single_df_ci['h']


manuscript_plot(single_df, dff, dff_ci, drug, ax2, ax4, False)


plt.tight_layout()
#plt.show()
plt.savefig("poor_ci_fit.pdf")
