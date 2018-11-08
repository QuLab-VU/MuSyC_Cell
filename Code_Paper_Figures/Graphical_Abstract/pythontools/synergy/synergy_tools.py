import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.stats import f as f_test
import pandas as pd
import os
import subprocess

#################################
# Helper functions
#################################
def hill_1D(d, E0, E1, h, ec50, logspace=False):
    """
    1D hill equation. Returns effect at given dose d.
    """
    if logspace:
        ec50 = np.power(10.,ec50)
        d = np.power(10.,d)
    return E1 + (E0-E1) / (1. + (d/ec50)**h)

def hill_1D_inv(E, E0, E1, h, ec50, logspace=False):
    """
    Inversion of 1D hill equation. Returns dose that achieves effect e.
    """
#    if (e < min(e0,e1)) or (e > max(e0,e1)): return np.nan
    if logspace: ec50 = np.power(10., ec50)
    d = np.power((E0-E1)/(E-E1)-1.,1./h)*ec50
    if logspace: return np.log10(d)
    else: return d


def hill_2D(d1, d2, E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r, logspace=False):
    """
    This is the 2D function used for MuSyC, assuming no cooperativity synergy
    """
    U, A1, A2, A12 = getUA_2D(d1, d2, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r, logspace=logspace)
    return U*E0 + A1*E1 + A2*E2 + A12*E3
    
def getUA_2D(d1, d2, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r, logspace=False):
    """
    This is the 2D function used for MuSyC, assuming no cooperativity synergy
    U, A1, A2, and A12 must be multiplied by appropriate effects to get surface height
    """
    if logspace:
        d1 = np.power(10.,d1)
        d2 = np.power(10.,d2)
    h1 = 1.*h1
    h2 = 1.*h2
    alpha1 = 1.*alpha1
    alpha2 = 1.*alpha2
    r1 = 1.*r1
    r1r = 1.*r1r
    r2 = 1.*r2
    r2r = 1.*r2r
    
    U = r1r*r2r*(r1*(alpha2*d1)**h1 + r1r + r2*(alpha1*d2)**h2 + r2r)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)

    A1 = r1*r2r*(d1**h1*r1*(alpha2*d1)**h1 + d1**h1*r1r + d1**h1*r2r + d2**h2*r2*(alpha2*d1)**h1)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)

    A2 = r1r*r2*(d1**h1*r1*(alpha1*d2)**h2 + d2**h2*r1r + d2**h2*r2*(alpha1*d2)**h2 + d2**h2*r2r)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)

    A12 = r1*r2*(d1**h1*r1*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r2r*(alpha1*d2)**h2 + d2**h2*r1r*(alpha2*d1)**h1 + d2**h2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1)/(d1**h1*r1**2*r2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d1**h1*r1**2*r2r*(alpha2*d1)**h1 + d1**h1*r1*r1r*r2*(alpha1*d2)**h2 + d1**h1*r1*r1r*r2r + d1**h1*r1*r2*r2r*(alpha1*d2)**h2 + d1**h1*r1*r2r**2 + d2**h2*r1*r1r*r2*(alpha2*d1)**h1 + d2**h2*r1*r2**2*(alpha1*d2)**h2*(alpha2*d1)**h1 + d2**h2*r1*r2*r2r*(alpha2*d1)**h1 + d2**h2*r1r**2*r2 + d2**h2*r1r*r2**2*(alpha1*d2)**h2 + d2**h2*r1r*r2*r2r + r1*r1r*r2r*(alpha2*d1)**h1 + r1r**2*r2r + r1r*r2*r2r*(alpha1*d2)**h2 + r1r*r2r**2)
    
    ### Or I couldj ust do A12 = 1 - (U+A1+A2)
    
    return U, A1, A2, A12
    
def u_hill(d1, d2, E1, E2, h1, h2, C1, C2):
    """
    From "Theory of synergistic effects: Hill-type response surfaces as 'null-interaction' models for mixtures" - Michael Schindler
    
    E - u_hill = 0 : Additive
    E - u_hill > 0 : Synergistic
    E - u_hill < 0 : Antagonistic
    """
    m1 = (d1/C1)
    m2 = (d2/C2)
    
    y = (h1*m1 + h2*m2) / (m1+m2)
    u_max = (E1*m1 + E2*m2) / (m1+m2)
    power = np.power(m1+m2, y)
    
    return u_max * power / (1. + power)

def zimmer_effective_dose(d1, d2, C1, C2, h1, h2, a12, a21, logspace=False):
    """
    From Zimmer's effective dose model
    May return NaNs if a12 or a21 are too negative
    
    Returns the effective dose surface value
    """
    if logspace:
        d1 = np.power(10., d1)
        d2 = np.power(10., d2)
        C1 = np.power(10., C1)
        C2 = np.power(10., C2)
        
    A = d2 + C2*(a21+1) + d2*a12
    B = d2*C1 + C1*C2 + a12*d2*C1 - d1*(d2+C2*(a21+1))
    C = -d1*(d2*C1 + C1*C2)

    d1p = (-B + np.sqrt(np.power(B,2.) - 4*A*C)) / (2.*A)
    d2p = d2 / (1. + a21 / (1. + C1/d1p))
    
    # Check a12 and a21 to make sure they are in a nice defined range
    return hill_1D(d1p, 1., 0., h1, C1) * hill_1D(d2p, 1., 0., h2, C2)

def rates(ec50, h, logspace=False, r1=1.):
    """
    U -> A, at rate r1*d**h
    A -> U, at rate r1r
    The dose d at which U=A=1/2 (half maximal effect) is called ec50
    
    At equilibrium
    r1*ec50**h = r1r
    r1r/r1 = ec50**h
    Assumes r1 = 1.
    returns r1=1, r1r=ec50**h
    """
    if logspace: ec50 = np.power(10.,ec50)
    r1r = np.power(ec50, h)
    return r1, r1r
    
def get_ec50(r1, r1r, h, logspace=False):
    """
    Returns ec50
    """
    ec50 = np.power(r1r/r1, 1./h)
    if logspace: return np.log10(ec50)
    return ec50
    
def fit_hill(d, E, sigma=None, logspace=False, E0=None, E1=None, E0_bounds=(-np.inf,np.inf), E1_bounds=(-np.inf,np.inf), h_bounds=(0,np.inf), ec50_bounds=(-np.inf,np.inf)):
    """
    Fits a dose response curve with doses given by d, effect given by E.
    
    Curve fit is always calculated based on log-transformed dose. If d is
    linear (logspace=False), d will be transformed into logspace, but ec50
    will be returned in linear space again
    
    Default: bounds = ([-np.inf, -np.inf, 0, -np.inf], np.inf)

    returns E0, E1, h, ec50, as ell as 1 x std of estimate of each parameter
    """
    if not logspace: d = np.log10(d)
    
    if E0 is not None and E1 is not None:
        f = lambda d, h, ec50: hill_1D(d, E0, E1, h, ec50, logspace=True)
        popt1, pcov = curve_fit(f, d, E, sigma=sigma, bounds=zip(h_bounds, ec50_bounds))
        h, ec50 = popt1

    elif E0 is not None:
        f = lambda d, E1, h, ec50: hill_1D(d, E0, E1, h, ec50, logspace=True)
        popt1, pcov = curve_fit(f, d, E, sigma=sigma, bounds=zip(E1_bounds, h_bounds, ec50_bounds))
        E1, h, ec50 = popt1
            
    elif E1 is not None:
        f = lambda d, E0, h, ec50: hill_1D(d, E0, E1, h, ec50, logspace=True)
        popt1, pcov = curve_fit(f, d, E, sigma=sigma, bounds=zip(E0_bounds, h_bounds, ec50_bounds))
        E0, h, ec50 = popt1
        
    else:
        f = lambda d, E0, E1, h, ec50: hill_1D(d, E0, E1, h, ec50, logspace=True)
        popt1, pcov = curve_fit(f, d, E, sigma=sigma, bounds=zip(E0_bounds, E1_bounds, h_bounds, ec50_bounds))
        E0, E1, h, ec50 = popt1

    perr = np.sqrt(np.diag(pcov))
        
    # From Christian
    p=1 #pvalue
    if sum(np.isnan(popt1))==0:
        # DIP rate fit
        E_fit = hill_1D(d, E0, E1, h, ec50)
        # F test vs flat linear "no effect" fit
        ssq_model = ((E_fit - E) ** 2).sum()
        ssq_null = ((np.mean(E) - E) ** 2).sum()
        df = len(d) - 2
        if E0 is None: df -= 1
        if E1 is None: df -= 1
        f_ratio = (ssq_null-ssq_model)/(ssq_model/df)
        p = 1 - f_test.cdf(f_ratio, 1, df)

    if not logspace:
        ec50 = np.power(10., ec50)
    return E0, E1, h, ec50, perr, p
    
def fit_hill_ci(d, E):
    """
    Fits a dose response curve with assuming that effect scales from 1
    to 0, using the median-effect equation from Chou-Talalay.
    
    d must be in linear space
    
    Requires 0 < e < 1: returns None if this is not satisfied
    
    returns 1, 0, h, ec50 (in linear space)
    """
    E = E.copy()
    E[(E < 0) | (E > 1)]=np.nan
    fU = E
    fA = 1-E

    C2 = linregress(np.log(d),np.log(fA/fU))
    h = C2.slope
    ec50 = np.exp(-C2.intercept/h)
    
    return 1., 0., h, ec50
    
def get_params(df, swap_drugs=False, logspace=False):
    E0 = np.float(df['E0'])
    E1 = np.float(df['E1'])
    E2 = np.float(df['E2'])
    E3 = np.float(df['E3'])

    h1 = np.float(df['h1'])
    h2 = np.float(df['h2'])

    r1 = np.float(df['r1'])
    r2 = np.float(df['r2'])
    
    C1 = np.float(df['log_C1'])
    C2 = np.float(df['log_C2'])
    
    d1min = np.float(df['min_conc_d1'])
    d2min = np.float(df['min_conc_d2'])
    
    d1max = np.float(df['max_conc_d1'])
    d2max = np.float(df['max_conc_d2'])
    
    if logspace:
        d1min = np.log10(d1min)
        d2min = np.log10(d2min)
        d1max = np.log10(d1max)
        d2max = np.log10(d2max)
        r1r = r1*np.power(np.power(10.,C1),h1)
        r2r = r2*np.power(np.power(10.,C2),h2)
    else:
        C1 = np.power(10.,C1)
        C2 = np.power(10.,C2)
        r1r = r1*np.power(C1,h1)
        r2r = r2*np.power(C2,h2)

    alpha1 = np.power(10.,np.float(df['log_alpha1']))
    alpha2 = np.power(10.,np.float(df['log_alpha2']))
    if (np.mean((E1,E2)) < E0):
        beta = np.float(df['beta'])/(E0-np.min((E1,E2)))
        beta_obs = np.float(df['beta_obs'])/(E0-df[['E1_obs','E2_obs']].min())
    else:
        beta = np.float(df['beta'])/(np.max((E1,E2)-E0))
        beta_obs = np.float(df['beta_obs'])/(df[['E1_obs','E2_obs']].max()-E0)
    

    
    model_level = str(df['model_level'])
    
    drug1_name = str(df['drug1_name'])
    drug2_name = str(df['drug2_name'])

    expt = str(df['expt'])
    
    if swap_drugs:
        atemp = alpha1
        ctemp = C1
        rtemp = r1
        r1rtemp = r1r
        e1temp = E1
        htemp = h1
        d1min_temp = d1min
        d1max_temp = d1max
        
        alpha1 = alpha2
        C1 = C2
        r1 = r2
        r1r = r2r
        E1 = E2
        h1 = h2
        d1min = d2min
        d1max = d2max
        
        alpha2 = atemp
        C2 = ctemp
        r2 = r1temp
        r2r = r1rtemp
        E2 = e1temp
        h2 = htemp
        d2min = d1min_temp
        d2max = d1max_temp
        
        drug2_name_temp = drug2_name
        drug2_name = drug1_name
        drug1_name = drug2_name_temp
    
    return E0, E1, E2, E3, C1, C2, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r, beta, beta_obs, d1min, d1max, d2min, d2max, model_level, drug1_name, drug2_name, expt

def get_example_params(E0=1., E1=0., E2=0., E3=0., alpha1=1.,      \
                       alpha2=1., d1min=-8., d2min=-8., d1max=-3., \
                       d2max=-3., C1=None, C2=None, h1=1., h2=1.,  \
                       r1=1., r2=1.):
    """
    Generates example parameters.
    Returns E0, E1, E2, E3, alpha1, alpha2, d1min, d1max, d2min, d2max, C1, C2, h1, h2, r1, r1r, r2, r2r
    """
    if C1 is None: C1 = (d1min + d1max)/2.
    if C2 is None: C2 = (d2min + d2max)/2.
    r1,r1r = rates(C1, h1, logspace=True, r1=r1)
    r2,r2r = rates(C2, h2, logspace=True, r1=r2)
    return E0, E1, E2, E3, alpha1, alpha2, d1min, d1max, d2min, d2max, C1, C2, h1, h2, r1, r1r, r2, r2r
    
def get_test_scatter_points(DD1, DD2, E, ndivs=1, reps=3, logspace=False):
    if logspace:
        DD1 = np.power(10.,DD1)
        DD2 = np.power(10.,DD2)

    nd2, nd1 = DD1.shape

    l = np.prod(DD1.shape)
    
    d1 = DD1.reshape((l,))
    d2 = DD2.reshape((l,))
    E = E.reshape((l,))
    
    emax = np.max(E)
    emin = np.min(E)
    erange = np.abs(emax-emin)

    df = pd.DataFrame(index = ["drug1.conc","drug2.conc","effect","effect.95ci"])
    
    d1i = -1
    d2i = 0
    for i in range(l):
        d1i += 1
        if d1i >= nd1:
            d1i=0
            d2i += 1
            if d2i > nd2: d2i=0
        
        if (d1i % ndivs != 0) or (d2i % ndivs != 0): continue
        for j in range(reps):
            r = (np.random.rand()*2-1)*erange/10.
            df[df.shape[1]] = [d1[i], d2[i], E[i]+r, np.random.rand()*erange/5.]
    
    return df.transpose()

#################################
# Synergy metrics
#################################
def combination_index(d1,d2,E):
    #Remove the values > 1
    #pv_sub.reset_index()
    pv_sub = pd.DataFrame([])
    pv_sub['per_via']=E
    pv_sub['drug1.conc']=d1
    pv_sub['drug2.conc']=d2
    #Remove the values > 1
    pv_sub = pv_sub[pv_sub['per_via']<1];
    #pv_sub.reset_index()
    d1_single = pv_sub[pv_sub['drug2.conc']==0]
    d2_single = pv_sub[pv_sub['drug1.conc']==0]
    d1_t = d1_single['drug1.conc'].values
    d2_t = d2_single['drug2.conc'].values
    fA1 = 1- d1_single['per_via'].values
    fA2 = 1- d2_single['per_via'].values
    fA1 = fA1[d1_t!=0]
    d1_t = d1_t[d1_t!=0]
    fA2 = fA2[d2_t!=0]
    d2_t = d2_t[d2_t!=0]
    
    #Fit the data to a line
    C1 = linregress(np.log(d1_t),np.log(fA1/(1-fA1)))
    C2 = linregress(np.log(d2_t),np.log(fA2/(1-fA2)))
    m1 = C1.slope
    m2 = C2.slope
    Dm1 = np.exp(C1.intercept/(-m1))
    Dm2 = np.exp(C2.intercept/(-m2))
    
    #Calculate the CI
    def d1_ci(fA):
        return (fA/(1-fA))**(1/m1)*Dm1
    def d2_ci(fA):
        return (fA/(1-fA))**(1/m2)*Dm2
    
    #The effect of the combination and the respective drug concentrations
    combin_effect = 1-pv_sub['per_via'][np.logical_and(pv_sub['drug1.conc']!=0 ,pv_sub['drug2.conc']!=0)].values
    d1_full = pv_sub['drug1.conc'][np.logical_and(pv_sub['drug1.conc']!=0 ,pv_sub['drug2.conc']!=0)].values
    d2_full = pv_sub['drug2.conc'][np.logical_and(pv_sub['drug1.conc']!=0 ,pv_sub['drug2.conc']!=0)].values
    CI = np.zeros((len(combin_effect),1));
    for e,eff in enumerate(combin_effect):
        CI[e] = d1_full[e]/d1_ci(eff) + d2_full[e]/d2_ci(eff)
    return CI

def combination_index_fit(d1,d2,E):
    #Remove the values > 1
    #pv_sub.reset_index()
    pv_sub = pd.DataFrame([])
    pv_sub['per_via']=E
    pv_sub['drug1.conc']=d1
    pv_sub['drug2.conc']=d2
    #Remove the values > 1
    pv_sub = pv_sub[pv_sub['per_via']<1];
    #pv_sub.reset_index()
    d1_single = pv_sub[pv_sub['drug2.conc']==0]
    d2_single = pv_sub[pv_sub['drug1.conc']==0]
    d1_t = d1_single['drug1.conc'].values
    d2_t = d2_single['drug2.conc'].values
    fA1 = 1- d1_single['per_via'].values
    fA2 = 1- d2_single['per_via'].values
    fA1 = fA1[d1_t!=0]
    d1_t = d1_t[d1_t!=0]
    fA2 = fA2[d2_t!=0]
    d2_t = d2_t[d2_t!=0]
    
    #Fit the data to a line
    C1 = linregress(np.log(d1_t),np.log(fA1/(1-fA1)))
    C2 = linregress(np.log(d2_t),np.log(fA2/(1-fA2)))
    m1 = C1.slope
    m2 = C2.slope
    Dm1 = np.exp(C1.intercept/(-m1))
    Dm2 = np.exp(C2.intercept/(-m2))
    
    #Calculate the CI
    def d1_ci(fA):
        return (fA/(1-fA))**(1/m1)*Dm1
    def d2_ci(fA):
        return (fA/(1-fA))**(1/m2)*Dm2
    
    #The effect of the combination and the respective drug concentrations
    combin_effect = 1-pv_sub['per_via'][np.logical_and(pv_sub['drug1.conc']!=0 ,pv_sub['drug2.conc']!=0)].values
    d1_full = pv_sub['drug1.conc'][np.logical_and(pv_sub['drug1.conc']!=0 ,pv_sub['drug2.conc']!=0)].values
    d2_full = pv_sub['drug2.conc'][np.logical_and(pv_sub['drug1.conc']!=0 ,pv_sub['drug2.conc']!=0)].values
    CI = np.zeros((len(combin_effect),1));
    for e,eff in enumerate(combin_effect):
        CI[e] = d1_full[e]/d1_ci(eff) + d2_full[e]/d2_ci(eff)
        
    return Dm1,Dm2,m1,m2


def loewe(d1, d2, E, E0, E1, E2, h1, h2, C1, C2, logspace=False):
    """
    Loewe = 1 : Additive
    Loewe > 1 : Antagonistic
    Loewe < 1 : Synergistic

    To get CI, use h and ec50 from fit_hill_ci, and use e0=1 and emax=0
    """
    if logspace:
        d1 = np.power(10., d1)
        d2 = np.power(10., d2)
        C1 = np.power(10., C1)
        C2 = np.power(10., C2)
    
    D1C = hill_1D_inv(E, E0, E1, h1, C1) # concentration of drug 1 alone which gets effect E
    D2C = hill_1D_inv(E, E0, E2, h2, C2) # concentration of drug 2 alone which gets effect E
    return d1/D1C + d2/D2C

def bliss(e_drug1_alone, e_drug2_alone, e_combination):
    """
    Bliss = 0 : Additive
    Bliss > 0 : Synergistic
    Bliss < 0 : Antagonistic 
    
    Unreliable behavior unless 0 <= e <= 1 for all arguments.
    """
    return e_drug1_alone*e_drug2_alone - e_combination
    
def hsa(e_drug1_alone, e_drug2_alone, e_combination, decreasing=True):
    """
    HSA = 0 : Additive
    HSA > 0 : Synergistic
    HSA < 0 : Antagonistic 
    
    If decreasing=True, assumes that the stronger effect is a lower value
    If decreasing=False, assumes that the stronger effect is a higher value
    """
    if len(e_drug1_alone)==1:
        if decreasing: return min(e_drug1_alone, e_drug2_alone) - e_combination
        else: return e_combination - max(e_drug1_alone, e_drug2_alone)
    else: 
        DD1 =e_drug1_alone
        DD2 =e_drug2_alone
        if decreasing:
            emx = 0
            for i in range(len(DD1)):
                emx = emx + min((e_combination[(DD1==DD1[i])&(DD2==0)],e_combination[(DD1==0)&(DD2==DD2[i])]))-e_combination[i]
            return emx
        else:
            emx = 0
            for i in range(len(DD1)):
                emx = emx - max((e_combination[(DD1==DD1[i])&(DD2==0)],e_combination[(DD1==0)&(DD2==DD2[i])]))+e_combination[i]
            return emx        
def schindler(d1, d2, E, E0, E1, E2, h1, h2, C1, C2, logspace=False):
    """
    From "Theory of synergistic effects: Hill-type response surfaces as 'null-interaction' models for mixtures" - Michael Schindler
    
    schin = 0 : Additive
    schin > 0 : Synergistic
    schin < 0 : Antagonistic
    """
    if logspace:
        d1 = np.power(10., d1)
        d2 = np.power(10., d2)
        C1 = np.power(10., C1)
        C2 = np.power(10., C2)
    return (E0-E) - u_hill(d1, d2, E0-E1, E0-E2, h1, h2, C1, C2)

def zimmer(DD1, DD2, E, logspace=False, sigma=None, vector_input=False):
    """
    Returns a12, a21, C1, C2, h1, and h2
    """
    if logspace:
        DD1 = np.power(DD1,10.)
        DD2 = np.power(DD2,10.)
    
    if not vector_input:
        d1 = DD1.reshape((-1,))
        d2 = DD2.reshape((-1,))
        E = E.reshape((-1,))
        
    else:
        d1 = DD1
        d2 = DD2
        E = E

    xdata = np.vstack((d1,d2))
        
    f = lambda d, C1, C2, h1, h2, a12, a21: zimmer_effective_dose(d[0], d[1], C1, C2, h1, h2, a12, a21, logspace=False)
    popt1, pcov = curve_fit(f, xdata, E, sigma=sigma, p0=[np.power(10.,-1),np.power(10.,-1),1,1,1,1],maxfev=10**4)
    C1, C2, h1, h2, a12, a21 = popt1
    perr = np.sqrt(np.diag(pcov))
    return C1, C2, h1, h2, a12, a21, perr    

def zimmer_old(DD1, DD2, E, C1, C2, h1, h2, logspace=False, sigma=None,vector_input=False):
    """
    Returns a12, a21 (requires C1, C2, h1, and h2)
    """

    if vector_input:
        d1 = DD1.reshape((np.prod(DD1.shape),))
        d2 = DD2.reshape((np.prod(DD2.shape),))
        E = E.reshape((np.prod(E.shape),))
    else:
        d1 = DD1
        d2 = DD2
        E = E
        
    if logspace:
        d1 = np.power(10., d1)
        d2 = np.power(10., d2)
        C1 = np.power(10., C1)
        C2 = np.power(10., C2)
        
    xdata = np.vstack((d1,d2))
    
    f = lambda d, a12, a21: zimmer_effective_dose(d[0], d[1], C1, C2, h1, h2, a12, a21, logspace=False)
    try:
        popt1, pcov = curve_fit(f, xdata, E, sigma=sigma,p0=[1,1],maxfev=10**6)
        perr = np.sqrt(np.diag(pcov))
    except RuntimeError:
        popt1 = [np.nan,np.nan]
        perr = [np.nan,np.nan]
    a12, a21 = popt1
    return a12, a21, perr


def zip(DD1,DD2,E,logspace=False,drug1_name='drug1',drug2_name='drug2',vector_input=False):
    if logspace:
        DD1 = np.power(10.,DD1)
        DD2 = np.power(10.,DD2)
    if not os.path.isdir('zip'):
        os.mkdir('zip')
        
    if vector_input:
        x1 = np.sort(np.unique(DD1))
        x2 = np.sort(np.unique(DD2))   
        xx1,xx2 = np.meshgrid(x1,x2)
        E_tmp = np.zeros((len(x1),len(x2)))
        for e1,i in enumerate(x1):
            for e2,j in enumerate(x2):
                E_tmp[e1,e2] = np.mean(E[(DD1==i)&(DD2==j)])
                
        E = E_tmp
    #Need three tables to pass to script in Delta_calculation.R
    #Plate_mat is the treated matrix
    #conc_range is the doses for the drugs
    #pair.list is meta data for the drugs
    #See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4759128/ for more info and examples
    plate_mat  = 1-pd.DataFrame(E).T
    plate_mat.to_csv('zip/plate_mat.csv',index=False)

    conc_range = pd.DataFrame([])
    if vector_input:
        conc_range['drug1name']= np.sort(np.unique(DD1))
        conc_range['drug2name']= np.sort(np.unique(DD2))    
    else:    
        conc_range['drug1name'] = DD1[0,:]
        conc_range['drug2name'] = DD2[:,0]
    conc_range.to_csv('zip/conc_range.csv',index=False)    

    pd.DataFrame([{'index':1,'drug1':drug1_name,'drug2':drug2_name,'cell.line':'cellline','plate':1}]).to_csv('zip/pair_list.csv')
        
    subprocess.call ("zip/calcZip.R")
    if os.path.exists('zip/delta_score.csv'):
        t = pd.read_csv('zip/delta_score.csv')
        os.remove('zip/delta_score.csv')
        return t['x'].values[0]
    else:
        return np.nan
   
   

def braid(d1,d2,e,logspace=False):
    if logspace:
        d1 = np.power(d1,10.)
        d2 = np.power(d2,10.)
    if not os.path.isdir('braid'):
        os.mkdir('braid')
    #https://www.nature.com/articles/srep25523
    #write effect in a dataframe
    df = pd.DataFrame()
    df['eff']=e
    df['d1.conc']=d1
    df['d2.conc']=d2
    
    df.to_csv('braid/effect.csv',index=False)
         
    subprocess.call ("braid/calcBraid.R")
    if os.path.exists('braid/calc_kappa.csv'):
        t = pd.read_csv('braid/calc_kappa.csv')
        os.remove('braid/calc_kappa.csv')
        return t['x'].values[0]
    else:
        return np.nan
   
def get_data(expt,drug1_name,drug2_name,target_cell_line):    
    #Read in the data into a pandas data frame
    data = pd.read_table(expt, delimiter=',')
    #Subset by target cell line
    data['drug1']=data['drug1'].str.lower()
    data['drug2']=data['drug2'].str.lower()
    data['sample']=data['sample'].str.upper()

    sub_data = data[data['sample']==target_cell_line]
    #data = data.reset_index()
    #control data for drug 1.  Either
    indx    = ((sub_data['drug1']==drug1_name) & (sub_data['drug2.conc']==0))  & (sub_data['drug1.conc']!=0)
    d1      = sub_data[indx]['drug1.conc'].values
    d2      = sub_data[indx]['drug2.conc'].values
    dip     = sub_data[indx]['effect'].values
    dip_95ci= sub_data[indx]['effect.95ci'].values
    
    indx    = ((sub_data['drug2']==drug1_name) & (sub_data['drug1.conc']==0))  & (sub_data['drug2.conc']!=0)
    d1      = np.concatenate([d1,sub_data[indx]['drug2.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug1.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    #control for drug 2
    indx    = ((sub_data['drug1']==drug2_name) & (sub_data['drug2.conc']==0))  & (sub_data['drug1.conc']!=0)
    d1      = np.concatenate([d1,sub_data[indx]['drug2.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug1.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    indx    = ((sub_data['drug2']==drug2_name) & (sub_data['drug1.conc']==0))  & (sub_data['drug2.conc']!=0)
    d1      = np.concatenate([d1,sub_data[indx]['drug1.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug2.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    #Combination experiment
    indx    = ((sub_data['drug1']==drug1_name) & (sub_data['drug2']==drug2_name)) & ((sub_data['drug1.conc']!=0) & (sub_data['drug2.conc']!=0)) 
    d1      = np.concatenate([d1,sub_data[indx]['drug1.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug2.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    indx    = ((sub_data['drug2']==drug1_name) & (sub_data['drug1']==drug2_name)) & ((sub_data['drug1.conc']!=0) & (sub_data['drug2.conc']!=0)) 
    d1      = np.concatenate([d1,sub_data[indx]['drug2.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug1.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    #The double control condition
    indx    = ((sub_data['drug1.conc']==0) & (sub_data['drug2.conc']==0))
    d1      = np.concatenate([d1,sub_data[indx]['drug1.conc'].values])
    d2      = np.concatenate([d2,sub_data[indx]['drug2.conc'].values])
    dip     = np.concatenate([dip,sub_data[indx]['effect'].values])
    dip_95ci= np.concatenate([dip_95ci,sub_data[indx]['effect.95ci'].values])
    
    #Set as standard deviation
    dip_sd = dip_95ci/(2*1.96)
    
    #Remove nan values
    d1      = d1[~np.isnan(dip)]
    d2      = d2[~np.isnan(dip)]
    dip_sd  = dip_sd[~np.isnan(dip)]
    dip     = dip[~np.isnan(dip)]
    
    return d1,d2,dip,dip_sd

    
if __name__ == "__main__":
    E0, E1, E2, E3, alpha1, alpha2, d1min, d1max, d2min, d2max, C1, C2, h1, h2, r1, r1r, r2, r2r = get_example_params(E2=0.5)
    d1 = np.logspace(d1min,d1max,10)
    d2 = np.logspace(d2min,d2max,10)
    DD1, DD2 = np.meshgrid(d1,d2)
    
    c1 = np.power(10.,C1)
    c2 = np.power(10.,C2)
    
    E = hill_2D(DD1, DD2, E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
    l = loewe(DD1, DD2, E, E0, E1, E2, h1, h2, c1, c2)
    
    E1_alone = hill_1D(DD1, E0, E1, h1, c1)
    E2_alone = hill_1D(DD2, E0, E2, h2, c2)
    b = bliss(E1_alone, E2_alone, E)
    
    s = schindler(DD1, DD2, E, E0, E1, E2, h1, h2, c1, c2)
    
    import synergy_plots
    import pandas as pd
    
    n_points = DD1.shape[0]*DD2.shape[1]
    df = pd.DataFrame(index = range(n_points))
    df['drug1.conc'] = DD1.reshape((n_points,))
    df['drug2.conc'] = DD2.reshape((n_points,))
    df['effect'] = E.reshape((n_points,))
    df['effect.95ci'] = np.nan
    for i in df.index:
        df.loc[i,'effect'] = df.loc[i,'effect']-(np.random.rand()-0.5)/4.
        df.loc[i,'effect.95ci'] =(np.random.rand()-0.5)/2.
    
    synergy_plots.plot_surface(d1min, d1max, d2min, d2max, E0, E1, E2, E3, h1, h2, r1, r1r, r2, r2r, alpha1, alpha2, d1buffer=3, d2buffer=2, scatter_points=df, zslice=0.2)


