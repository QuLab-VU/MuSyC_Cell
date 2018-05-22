import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.stats import linregress, pearsonr, spearmanr
from scipy.stats import kruskal

def mag_max(x):
    x = np.asarray(x) # Make sure x is a numpy array
    i = np.argmax(np.abs(x), axis=1) # Get the index of the max absolute value of every row
    return np.diag(x[:,i])

def quadrant(log_alpha1, log_alpha2, beta):
    if log_alpha1 > 0: x = 1
#    else: x = 5
    else: x = 1
    
    if beta > 0:
        if log_alpha2 > 0: return x
        else: return 1+x
    else:
        if log_alpha2 > 0: return 3+x
        else: return 2+x
        
def quadrant_name(x):
    if x==1: return "Q.I"
    if x==2: return "Q.II"
    if x==3: return "Q.III"
    if x==4: return "Q.IV"
    return "?err"



master_drug_results = pd.read_csv("../../Data/MasterResults_PC9_filtered.csv")
syn_params_dict = dict()
for col in ['drug','bliss','loewe','ci']: syn_params_dict[col] = []

drug_params_dict = dict()
for col in ['drug_name','log_alpha_1','log_alpha_2','beta','hill','hill_dip','ci_hill','log_ec50','emax','dip_edrug','max_conc']: drug_params_dict[col] = []

blc18 = pd.read_csv("OtherMetrics/results/HTS018_bliss_loewe_ci.csv")
blc15 = pd.read_csv("OtherMetrics/results/HTS015_bliss_loewe_ci.csv")
blccella = pd.read_csv("OtherMetrics/results/cellavista_bliss_loewe_ci.csv")
lodip18 = pd.read_csv("loewe_from_dip/HTS018_loewe.csv")
lodip15 = pd.read_csv("loewe_from_dip/HTS015_loewe.csv")
lodipcella = pd.read_csv("loewe_from_dip/cellavista_loewe.csv")



hill_params_18 = pd.read_csv("OtherMetrics/single_responses/18/parameters.csv", index_col=0)
hill_params_15 = pd.read_csv("OtherMetrics/single_responses/15/parameters.csv", index_col=0)
hill_params_cella = pd.read_csv("OtherMetrics/single_responses/cellavista/parameters.csv", index_col=0)
hill_params_18.index = [i.lower() for i in hill_params_18.index]
hill_params_15.index = [i.lower() for i in hill_params_15.index]
hill_params_cella.index = [i.lower() for i in hill_params_cella.index]
for i in master_drug_results.index:
    if master_drug_results.loc[i,'expt'] == "HTS018_rates.csv":
        drug_id = "drug2"
        osi_id = "drug1"
        drug = master_drug_results.loc[i,"%s_name"%drug_id]
        df = blc18.loc[(blc18[drug_id]==drug) & (blc18['drug1.conc'] > 0) & (blc18['drug2.conc'] > 0)]
        df_loewe = lodip18.loc[(lodip18[drug_id]==drug) & (lodip18['drug1.conc'] > 0) & (lodip18['drug2.conc'] > 0)]
        drug_params_dict['drug_name'].append(drug)
        drug_params_dict['log_alpha_1'].append(master_drug_results.loc[i,"log_alpha1"])
        drug_params_dict['log_alpha_2'].append(master_drug_results.loc[i,"log_alpha2"])
        drug_params_dict['hill_dip'].append(master_drug_results.loc[i,"h2"])
        drug_params_dict['beta'].append(master_drug_results.loc[i,"beta"])
        drug_params_dict['log_ec50'].append(master_drug_results.loc[i,"log_C2"])
        drug_params_dict['dip_edrug'].append(master_drug_results.loc[i,"E2"])
        drug_params_dict['max_conc'].append(np.log10(master_drug_results.loc[i,"max_conc_d2"]))
        drug_params_dict['hill'].append(hill_params_18.loc[drug,"h"])
        drug_params_dict['emax'].append(hill_params_18.loc[drug,"emax"])
        drug_params_dict['ci_hill'].append(hill_params_18.loc[drug,"h_CI"])

    elif master_drug_results.loc[i,'expt'] == "HTS015_017_Combined.csv":
        drug_id = "drug1"
        osi_id = "drug2"
        drug = master_drug_results.loc[i,"%s_name"%drug_id]
        df = blc15.loc[(blc15[drug_id]==drug) & (blc15['drug1.conc'] > 0) & (blc15['drug2.conc'] > 0)]
        df_loewe = lodip15.loc[(lodip15[drug_id]==drug) & (lodip15['drug1.conc'] > 0) & (lodip15['drug2.conc'] > 0)]
        drug_params_dict['drug_name'].append(drug)
        drug_params_dict['log_alpha_1'].append(master_drug_results.loc[i,"log_alpha2"])
        drug_params_dict['log_alpha_2'].append(master_drug_results.loc[i,"log_alpha1"])
        drug_params_dict['log_ec50'].append(master_drug_results.loc[i,"log_C1"])
        drug_params_dict['beta'].append(master_drug_results.loc[i,"beta"])
        drug_params_dict['hill_dip'].append(master_drug_results.loc[i,"h1"])
        drug_params_dict['ci_hill'].append(hill_params_15.loc[drug,"h_CI"])
        drug_params_dict['hill'].append(hill_params_15.loc[drug,"h"])
        drug_params_dict['emax'].append(hill_params_15.loc[drug,"emax"])
        drug_params_dict['dip_edrug'].append(master_drug_results.loc[i,"E1"])
        drug_params_dict['max_conc'].append(np.log10(master_drug_results.loc[i,"max_conc_d1"]))
    else:
        drug_id = "drug1"
        osi_id = "drug2"
        drug = master_drug_results.loc[i,"%s_name"%drug_id]
        df = blccella.loc[(blccella[drug_id]==drug) & (blccella['drug1.conc'] > 0) & (blccella['drug2.conc'] > 0)]
        df_loewe = lodipcella.loc[(lodipcella[drug_id]==drug) & (lodipcella['drug1.conc'] > 0) & (lodipcella['drug2.conc'] > 0)]
        drug_params_dict['drug_name'].append(drug)
        drug_params_dict['log_alpha_1'].append(master_drug_results.loc[i,"log_alpha2"])
        drug_params_dict['log_alpha_2'].append(master_drug_results.loc[i,"log_alpha1"])
        drug_params_dict['log_ec50'].append(master_drug_results.loc[i,"log_C1"])
        drug_params_dict['beta'].append(master_drug_results.loc[i,"beta"])
        drug_params_dict['ci_hill'].append(hill_params_cella.loc[drug,"h_CI"])
        drug_params_dict['hill_dip'].append(master_drug_results.loc[i,"h1"])
        drug_params_dict['hill'].append(hill_params_cella.loc[drug,"h"])
        drug_params_dict['emax'].append(hill_params_cella.loc[drug,"emax"])
        drug_params_dict['dip_edrug'].append(master_drug_results.loc[i,"E1"])
        drug_params_dict['max_conc'].append(np.log10(master_drug_results.loc[i,"max_conc_d1"]))
        
    for j in df.index:
        syn_params_dict['drug'].append(drug.lower())
        syn_params_dict['bliss'].append(df.loc[j,'bliss'])
#        if pd.isnull(df.loc[j, "loewe"]): syn_params_dict['loewe'].append(np.nan)
#        else: syn_params_dict['loewe'].append(-np.log10(df.loc[j, "loewe"]))
        if pd.isnull(df.loc[j, "CI"]): syn_params_dict['ci'].append(np.nan)
        else: syn_params_dict['ci'].append(-np.log10(df.loc[j, "CI"]))
    for j in df_loewe.index:
        if pd.isnull(df_loewe.loc[j, "loewe"]): syn_params_dict['loewe'].append(np.nan)
        else: syn_params_dict['loewe'].append(-np.log10(df_loewe.loc[j, "loewe"]))
    if drug=="dasatinib": syn_params_dict['loewe'].append(np.nan) # Somewhere we lost a single value of dasatinib from the DIP rate data
#    print drug, df.shape[0], df_loewe.shape[0]

syn_params = pd.DataFrame(syn_params_dict)
drug_params = pd.DataFrame(drug_params_dict)
drug_params.index = drug_params['drug_name']

bliss = syn_params.loc[~(pd.isnull(syn_params['bliss'])), ['drug','bliss']]
loewe = syn_params.loc[~(pd.isnull(syn_params['loewe'])), ['drug','loewe']]
chou = syn_params.loc[~(pd.isnull(syn_params['ci'])), ['drug','ci']]

#raise(ValueError)

quadrants = []
for i in drug_params.index:
    quadrants.append(quadrant(drug_params.loc[i,"log_alpha_1"], drug_params.loc[i,"log_alpha_2"], drug_params.loc[i,"beta"]))
drug_params['quadrant']=quadrants

drug_mechanisms = pd.read_csv("../../Data/Drugs_MOA_3-22-18.csv",index_col=0)
drug_params['class_order'] = [str(i).split('.')[0] for i in drug_mechanisms.loc[drug_params.index,"label"]]


drug_class_colors = dict()
drug_class_colors['3'] = "D62728" # Red
drug_class_colors['4'] = "9467BD" # Purple
drug_class_colors['2'] = "FF7F0E" # Orange
drug_class_colors['1'] = "1F77B4" # Blue
drug_class_colors['Other'] = "000000"

drug_params['color'] = [drug_class_colors[i] for i in drug_params['class_order']]
#drug_params['drug_name'] = drug_params.index

drug_params.sort_values(by=["quadrant", "class_order", "drug_name"], ascending=[True,True,True], inplace=True)

bliss_data = [bliss.loc[bliss['drug']==drug,'bliss'] for drug in drug_params.index]
loewe_data = [loewe.loc[loewe['drug']==drug,'loewe'] for drug in drug_params.index]
chou_data = list()
for drug in drug_params.index:
    if drug_params.loc[drug,'ci_hill']<0: chou_data.append([])
    else: chou_data.append(chou.loc[chou['drug']==drug,'ci'])

drug_params['bliss_median'] = np.nan
drug_params['loewe_median'] = np.nan
drug_params['ci_median'] = np.nan
drug_params['bliss_min'] = np.nan
drug_params['loewe_min'] = np.nan
drug_params['ci_min'] = np.nan
drug_params['bliss_max'] = np.nan
drug_params['loewe_max'] = np.nan
drug_params['ci_max'] = np.nan
for drug, bd, ld, cd in zip(drug_params.index, bliss_data, loewe_data,chou_data):
    drug_params.loc[drug,'bliss_median'] = np.median(bd)
    drug_params.loc[drug,'loewe_median'] = np.median(ld)
    drug_params.loc[drug,'ci_median'] = np.median(cd)
    drug_params.loc[drug,'bliss_min'] = np.min(bd)
    drug_params.loc[drug,'loewe_min'] = np.min(ld)
    try: drug_params.loc[drug,'ci_min'] = np.min(cd)
    except ValueError: drug_params.loc[drug,'ci_min']=np.nan
    drug_params.loc[drug,'bliss_max'] = np.max(bd)
    drug_params.loc[drug,'loewe_max'] = np.max(ld)
    try: drug_params.loc[drug,'ci_max'] = np.max(cd)
    except ValueError: drug_params.loc[drug,'ci_max']=np.nan
    

l = drug_params.shape[0]
drug_list = drug_params.index



drug_short_names = dict()
drug_short_names['beclomethasonedipropionate'] = 'beclo'
drug_short_names['bendroflumethiazide'] = 'bendro'
drug_short_names['cephalomannine'] = 'cephalo'
drug_short_names['tyrphostinag370'] = 'tyrphostin'


counts = drug_params['quadrant'].value_counts()
counts = counts[counts.index.sort_values()]


if True: ################################################## ******************** WE WANT THIS ONE *********************
    fig = plt.figure(figsize=(4,2))
    
    
    dp = drug_params
    
    x_ci = dp['ci_hill']
    y_ci = dp['ci_median']
    q_ci = dp['quadrant']
    x_ci = x_ci[~pd.isnull(y_ci)]
    y_ci = y_ci[~pd.isnull(y_ci)]
    
    x_lo = dp['hill_dip']
    y_lo = dp['loewe_median']
    q_lo = dp['quadrant']
    x_lo = x_lo[~pd.isnull(y_lo)]
    q_lo = q_lo[~pd.isnull(y_lo)]
    y_lo = y_lo[~pd.isnull(y_lo)]
    
    ax = fig.add_subplot(121)
    ax.scatter(x_lo, y_lo, alpha=0.8, s=20)
    ax.set_title("Loewe", fontsize=8)
    ax.set_xlabel("DIP h2", fontsize=8)
    ax.set_ylabel("Median Synergy", fontsize=8)
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    loewe_line = linregress(x_lo, y_lo)
    ax.plot(xlim, [loewe_line.intercept + xlim[0]*loewe_line.slope, loewe_line.intercept + xlim[1]*loewe_line.slope])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.tick_params(labelsize=8)
    
    
    ax = fig.add_subplot(122)
    ax.scatter(x_ci, y_ci, label="CI", alpha=0.8, s=20)
    ax.set_title("CI", fontsize=8)
    ax.set_xlabel("CI h2", fontsize=8)
    ylim = list(ax.get_ylim())
    ylim[1]=1
    xlim = ax.get_xlim()
    chou_line = linregress(x_ci, y_ci)
    ax.plot(xlim, [chou_line.intercept + xlim[0]*chou_line.slope, chou_line.intercept + xlim[1]*chou_line.slope])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.tick_params(labelsize=8)
    

    print "Pearson Loewe",pearsonr(x_lo,y_lo)[0],"CI",pearsonr(x_ci,y_ci)[0]
    print "Spearman Loewe",spearmanr(x_lo,y_lo)[0],"CI",spearmanr(x_ci,y_ci)[0]
    
    res = """
        Pvals:
        Pearson: Loewe=6.69297285854e-05,   CI=0.0228035585082
        Spearman Loewe=3.22853424366e-05,   CI=0.0245427460673
        
        cors:
        Pearson Loewe=-0.55850060934,        CI=-0.303858797903
        Spearman Loewe=-0.577865612648,      CI=-0.3002734108


    """
    

    plt.tight_layout()
    plt.savefig("Loewe_CI_Hill.pdf")



#raise(ValueError)
fig = plt.figure(figsize=(8,3))
gs = gridspec.GridSpec(5,counts.shape[0],height_ratios=[0.16,0.5,1,1,1], width_ratios = counts.values)

left = 0

loewe_quadrant = []
for i in range(counts.shape[0]): loewe_quadrant.append([])
ci_quadrant = []
for i in range(counts.shape[0]): ci_quadrant.append([])
bliss_quadrant = []
for i in range(counts.shape[0]): bliss_quadrant.append([])


for qi,q in enumerate(counts.index):



    # Drug class
    ax = plt.subplot(gs[0,qi])

    drug_labels = []
    colors = np.zeros((1,counts[q],3))
    for di,d in enumerate(drug_list[left:left+counts[q]]):
        if drug_short_names.has_key(d): drug_labels.append(drug_short_names[d])
        else: drug_labels.append(d)
        colors[0,di,:] = [ii / 255. for ii in bytearray.fromhex(drug_params.loc[d,"color"])]
        if di < counts[q]-1: plt.plot([di+0.5,di+0.5],[-0.5,0.5],'w-',lw=0.5,alpha=0.5)
        
#    print left, counts[q], len(drug_labels)
    real_x_ticks = [i for i in range(counts[q])]
    x_ticks = []
    x_tick_labels = []
    for t,d in zip(real_x_ticks, drug_labels):
        x_ticks.append(t-0.5)
        x_ticks.append(t)
        x_tick_labels.append(d)
        x_tick_labels.append("")
#        ax.set_xticks([i for i in range(counts[q])])
    ax.set_xticks(x_ticks)
    
    ax.set_xlim(-0.5,counts[q]-0.5)
    ax.set_ylim(-0.5,0.5)
    #ax.set_xticklabels(drug_labels, rotation=70, ha="left", fontsize=7)
    ax.set_xticklabels(x_tick_labels, rotation=65, ha="left", fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks([])
    ax.imshow(colors, aspect="auto")
    

    
    
    i = 0
    for tick in ax.xaxis.get_ticklines():
        if (i+1)%4 != 0: tick.set_markersize(0)
        else: tick.set_markersize(3.5)
        i += 1
    
    
    
    
    

    # Quadrant info
    ax = plt.subplot(gs[1,qi])

    ax.plot([0,0],[-1,1],'k-', lw=0.3)
    ax.plot([-1,1],[0,0],'k-', lw=0.3)
    if qi == 0: ax.fill_between([0,1],0,1,color="#FFCCCC")
    elif qi == 1: ax.fill_between([-1,0],0,1,color="#FFCCCC")
    elif qi == 2: ax.fill_between([-1,0],0,-1,color="#FFCCCC")
    elif qi == 3: ax.fill_between([0,1],0,-1,color="#FFCCCC")
    ax.set_ylim(-1.1,1.1)
    ax.set_xlim(-1.1,1.1)
    ax.set_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])



    # Loewe
    ax = plt.subplot(gs[2,qi])
    if q == 1: ax.set_ylabel("Loewe",fontsize=12)
    ax.boxplot(loewe_data[left:left+counts[q]],showfliers=False)
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.xaxis.tick_top()
    ax.yaxis.tick_right()
    xlim = ax.get_xlim()
#    ylim = ax.get_ylim()
    ylim = [-4, 4]
    if qi==3:
        ax.set_yticks([i/2. for i in ylim])
        ax.set_yticklabels(["Ant","Syn"])
    ax.fill_between([xlim[0],xlim[1]], ylim[1], 0, color='k', alpha=0.1)
    ax.set_ylim(ylim)
#    print "Loewe %s"%repr(ylim)
    for j in range(left, left+counts[q]): loewe_quadrant[qi].append(loewe_data[j].median())
    print qi,len(loewe_quadrant[qi])



    # CI
    ax = plt.subplot(gs[3,qi])
    if q ==1: ax.set_ylabel("CI",fontsize=12)
    ax.boxplot(chou_data[left:left+counts[q]],showfliers=False)
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.xaxis.tick_top()
    ax.yaxis.tick_right()
    xlim = ax.get_xlim()
    #ylim = list(ax.get_ylim())
    ylim = [-4, 4]
    #    ylim[1] = 0.5*ylim[1]
    #    ylim[0] = -1*ylim[1]
    #print "Chou %s"%repr(ylim)
    if qi==3:
        ax.set_yticks([i/2. for i in ylim])
        ax.set_yticklabels(["Ant","Syn"])
    ax.fill_between([xlim[0],xlim[1]], ylim[1], 0, color='k', alpha=0.1)
    ax.set_ylim(ylim)

    for j in range(left, left+counts[q]): ci_quadrant[qi].append(np.median(chou_data[j]))



    # Bliss
    ax = plt.subplot(gs[4,qi])
    if q==1: ax.set_ylabel("Bliss",fontsize=12)
    ax.boxplot(bliss_data[left:left+counts[q]],showfliers=False)
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.xaxis.tick_top()
    ax.yaxis.tick_right()
    xlim = ax.get_xlim()
#    ylim = ax.get_ylim()
    ylim = [-0.4, 0.3]
    if qi==3:
        ax.set_yticks([i/2. for i in ylim])
        ax.set_yticklabels(["Ant","Syn"])
    ax.fill_between([xlim[0],xlim[1]], ylim[1], 0, color='k', alpha=0.1)
    ax.set_ylim(ylim)
#    print "Bliss %s"%repr(ylim)

    for j in range(left, left+counts[q]): bliss_quadrant[qi].append(bliss_data[j].median())



    left += counts[q]
    
    loewe_quadrant[qi] = [i for i in loewe_quadrant[qi] if not pd.isnull(i)]
    ci_quadrant[qi] = [i for i in ci_quadrant[qi] if not pd.isnull(i)]
    bliss_quadrant[qi] = [i for i in bliss_quadrant[qi] if not pd.isnull(i)]


print "Kruskal-Wallis test on median Loewe across quadrants: %0.3e"%kruskal(*loewe_quadrant)[1]
print "Kruskal-Wallis test on median CI across quadrants: %0.3e"%kruskal(*ci_quadrant)[1]
print "Kruskal-Wallis test on median Bliss across quadrants: %0.3e"%kruskal(*bliss_quadrant)[1]


results_for_full_octants = """
Kruskal-Wallis test on median Loewe across quadrants: 3.672e-01      = 0.367
Kruskal-Wallis test on median CI across quadrants: 3.153e-01         = 0.315
Kruskal-Wallis test on median Bliss across quadrants: 4.578e-01      = 0.458
"""

results_for_quadrants = """
Kruskal-Wallis test on median Loewe across quadrants: 9.770e-01      = 0.977
Kruskal-Wallis test on median CI across quadrants: 3.153e-01         = 0.315 
Kruskal-Wallis test on median Bliss across quadrants: 3.183e-01      = 0.318
"""

plt.tight_layout()
plt.subplots_adjust(hspace=0.1, wspace=0.1)
plt.savefig("other_metrics_distributions.pdf")


