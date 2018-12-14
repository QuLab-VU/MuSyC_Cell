import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

def f(df, loewe_synergistic, loewe_antagonistic, loewe_und):
    for i in df.index:
        if df.loc[i,'drug1.conc'] == 0 or df.loc[i,'drug2.conc'] == 0: continue
        elif pd.isnull(df.loc[i,'loewe']): loewe_und.append(df.loc[i,'rate'])
        elif df.loc[i,'loewe'] < 1: loewe_synergistic.append(df.loc[i,'rate'])
        else: loewe_antagonistic.append(df.loc[i,'rate'])

blc18 = pd.read_csv("HTS018_loewe.csv")
blc15 = pd.read_csv("HTS015_loewe.csv")
blccella = pd.read_csv("cellavista_loewe.csv")


loewe_synergistic = []
loewe_antagonistic = []
loewe_und = []

f(blc18, loewe_synergistic, loewe_antagonistic, loewe_und)
f(blc15, loewe_synergistic, loewe_antagonistic, loewe_und)
f(blccella, loewe_synergistic, loewe_antagonistic, loewe_und)



fig = plt.figure(figsize=(2,3))
ax = fig.add_subplot(111)
ax.boxplot([loewe_synergistic,loewe_antagonistic, loewe_und], widths=0.8)
ax.set_ylabel(r"DIP rate ($s^{-1}$)", fontsize=8)
ax.set_xticklabels(['SYN','ANT','UND'], fontsize=8)
ax.tick_params(labelsize=8)
xlim = ax.get_xlim()
ax.plot(xlim,[0,0],'k--')
ax.set_xlim(xlim)
plt.tight_layout()
#plt.show()
plt.savefig("loewe_limitations.pdf")
