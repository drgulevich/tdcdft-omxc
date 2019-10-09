#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rc('font', size=14)
fig, ax = plt.subplots(2,1,figsize=(8,8))

data = np.loadtxt('out')

trange = data[:,0]
drange = data[:,1]

dt = trange[1]-trange[0]
tmax = trange[-1]

domega = 2.*np.pi/tmax
fd = dt*np.fft.fft(drange)
esd = (fd.conj()*fd).real # energy spectral density for a transient
omrange = domega*np.arange(len(esd))

ax[0].set_xlabel(r'$t$')
ax[0].set_ylabel(r'$d(t)$')
#ax[0].set_xlim([0,100.])
#ax[0].set_ylim([0,1.])
ax[0].plot(trange, drange, 'k')

omegacut = 3.0
indcut=int(round(omegacut/domega))
indmax = np.argmax(esd[1:indcut])+1

Eh = 11.271
print('# Warning! Check that you are using Eh=', Eh)
print('# omega at maximum: ',Eh*omrange[indmax])

print('Check Parseval theorem:')
print('np.sum(drange)*dt: ',np.sum(drange*drange)*dt)
print('np.sum(esd)*domega: ',np.sum(esd)*domega/(2.*np.pi))
print('2.*np.sum(esd[:indcut])*domega: ',2.*np.sum(esd[:indcut])*domega/(2.*np.pi))

ax[1].set_xlabel(r'$\omega,\rm\,meV$')
ax[1].set_xlim([0,30])
ax[1].set_ylabel(r'$|F[d(t)]|^2,\rm\,a.u.$')
ax[1].semilogy(Eh*omrange[1:indcut], esd[1:indcut], 'k')

plt.grid()
plt.tight_layout()
plt.show()
