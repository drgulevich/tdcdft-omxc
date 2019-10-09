#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rc('font', size=14)
fig, ax = plt.subplots(figsize=(9,6))
#ax.set_title("title")
ax.set_xlabel(r'$t$ (a.u.)')
ax.set_ylabel(r'$d(t)$ (a.u.)')
#ax.set_xlim([0,100.])
#ax.set_ylim([-2e-3,3e-3])


for f in sorted([f for f in os.listdir('.') if os.path.isfile(f)]):
    if f.startswith('out_a0') or f.startswith('out_m0'):
        print(f)
        data = np.loadtxt(f)
        ax.plot(data[:,0], data[:,1], alpha=0.4)
        #ax.scatter(data[:,0], data[:,1], s=2)


data = np.loadtxt('out')
ax.plot(data[:,0], data[:,1], 'k')

#ax.legend(loc='upper left')

plt.grid()
plt.tight_layout()
plt.show()
