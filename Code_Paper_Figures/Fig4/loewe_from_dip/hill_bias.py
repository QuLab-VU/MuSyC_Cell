import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

def getD2(d1, E0 = 1, E1 = 0, E2 = 0, E = 0.1, c1 = np.power(10.,-6.), c2 = np.power(10., -6.), h1=1., h2=1.):
    xxx = ((E - E0) + (E-E1)*(d1/c1)**h1) / (E2-E)
#    (d2/c2)**h2 = ((E - E0) + (E-E1)*(d1/c1)**h1) / (E2-E)
 #   h2 log(d2/c2) = xxx
 #   log(d2/c2) = xxx/h2
 #   d2/c2 = np.exp(log(xxx)/h2)
    
    return c2*np.exp(np.log(xxx)/h2)

def linear(y0, y1, x):
    x0 = min(x)
    x1 = max(x)
    m = (y1-y0)/(x1-x0)
    return y0 + m*(x-x0)

fig = plt.figure(figsize=(4,1.5))


# h=1
ax = fig.add_subplot(132)
x = np.linspace(0,getD2(0.),1000)
y = getD2(x)
ax.plot(x, y)
ax.fill_between(x, y, alpha=0.5)
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("h1 = h2 = 1", fontsize=8)
ax.set_xlabel("d1 conc. (linear)", fontsize=8)
ax.text(x[10],x[30],"Synergistic", fontsize=8)
ax.text(x[300],x[900],"Antagonistic", fontsize=8)
ax.set_xlim(0,max(x))
ax.set_ylim(0,max(x))



# h=2
h1 = 1.
h2 = 2.
ax = fig.add_subplot(133)
z = getD2(0.,h1=h1,h2=h2)
z2 = getD2(0.,h1=h2,h2=h1)
x = np.linspace(0,z2,1000)
x2 = np.linspace(0,z,1000)
ylinear = linear(z,0,x)
y = getD2(x, h1=h1, h2=h2)
ax.plot(x, y)
ax.plot([0,z2], [z,0],'k--')
ax.fill_between(x, y, alpha=0.5)
ax.fill_between(x, y, ylinear, alpha=0, hatch="/")
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("h1 = %0.1f, h2 = %0.1f"%(h1,h2), fontsize=8)
ax.set_xlabel("d1 conc. (linear)", fontsize=8)
ax.text(x[10],x2[30],"Synergistic", fontsize=8)
ax.text(x[300],x2[900],"Antagonistic", fontsize=8)
ax.set_xlim(0,z2)
ax.set_ylim(0,z)


# h=0.5
h1 = 1.
h2 = 0.5
ax = fig.add_subplot(131)
z = getD2(0.,h1=h1,h2=h2)
z2 = getD2(0.,h1=h2,h2=h1)
x = np.linspace(0,z2,1000)
x2 = np.linspace(0,z,1000)
ylinear = linear(z,0,x)
y = getD2(x, h1=h1, h2=h2)
ax.plot(x, y)
ax.plot([0,z2], [z,0],'k--')
ax.fill_between(x, y, alpha=0.5)
ax.fill_between(x, y, ylinear, alpha=0, hatch="/")
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("h1 = %0.1f, h2 = %0.1f"%(h1,h2), fontsize=8)
ax.set_xlabel("d1 conc. (linear)", fontsize=8)
ax.set_ylabel("d2 conc. (linear)", fontsize=8)
ax.text(x[10],x2[30],"Synergistic", fontsize=8)
ax.text(x[300],x2[900],"Antagonistic", fontsize=8)
ax.set_xlim(0,z2)
ax.set_ylim(0,z)

plt.tight_layout()

plt.savefig("hill_problem.pdf")
#plt.show()




