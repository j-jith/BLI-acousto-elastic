from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt
from scipy.optimize import curve_fit
import deal_vec
#import piecewise_polynomial_fit as ppf
#import vectfit as vf

# No. of data files
n_steps = 50

# Read data
print('Reading data ...')

kt_ratio = []
for k in range(n_steps):
    filename = 'ktr_op/ktr_{}.dat'.format(k)
    kt_ratio.append(deal_vec.read(filename))

kt_ratio = np.array(kt_ratio).T


# Compute correlation matrix 
print('Computing correlation matrix ...')
corr_mat = np.dot(kt_ratio.T, kt_ratio)

# Solve eigenvalue problem
print('Solving eigenvalue problem ...')
l, psi = np.linalg.eig(corr_mat)

# To check significance of eigenvalues
eig_sig = np.cumsum(l**2)/np.sum(l**2)

print('Computing POD modes and coefficients ...')
# Compute POD modes
phi = []
for k in range(len(psi)):
    phi.append(np.dot(kt_ratio, psi[:, k])) # POD modes
    phi[-1] = phi[-1]/np.linalg.norm(phi[-1])

phi = np.array(phi).T
# Compute POD coefficients
a = np.dot(phi.T, kt_ratio)

# Compute and plot error
print('Computing RMS error ...')
fig = plt.figure(1)
ax = fig.add_subplot(111)

for n_eig in range(40, 50):

    kt_recon = np.dot(phi[:, :n_eig],  a[:n_eig, :])

    err = np.array([])
    for k in range(n_steps):
        err = np.append(err, np.sqrt(np.mean((kt_ratio[:, k]-kt_recon[:, k])**2))/(np.max(kt_ratio[:, k]) - np.min(kt_ratio[:, k])))

    ax.semilogy(np.linspace(1, n_steps, n_steps), err, label='n = {}'.format(n_eig+1))

print('Done!')

ax.set_ylabel('Normalised RMS error')
ax.set_xlabel('Step')
ax.legend()
plt.show()


## Piecewise polynomial fit

## Fluid properties
rho = 286.13
c = 283.76
cp = 863.6
gamma = 1.28
kappa = 0.0754
mu = 5.546e-5

## Global frequency range
f_0 = 150 # Hz
f_1 = 750 # Hz
w_g = np.linspace(2*pi*f_0, 2*pi*f_1, n_steps) # radians

pieces = 2
w_g_split = np.array_split(w_g, pieces)
f_pieces = [ww[0]/2/pi for ww in w_g_split]
f_pieces.append(w_g[-1]/2/pi)

## Figures
fig_fit = plt.figure(1)
ax_fit = fig_fit.add_subplot(111)

fig_err = plt.figure(2)
ax_err = fig_err.add_subplot(111)

def func_fit(x, e0, e1, e2):
    return e0 + e1*x + e2/x

## Fitting loop
for ww in w_g_split:
    
    k02 = ww**2/c**2
    kv2 = -1j*ww*rho/mu
    kh2 = -1j*ww*rho*cp/kappa

    ii = np.argmin(np.abs(w_g - ww[0]))
    ii_1 = np.argmin(np.abs(w_g - ww[-1]))

    ff = a[0, ii:ii_1]
    ff = ff*np.sqrt(k02/kv2) + (gamma-1)*np.sqrt(k02/kh2)
    #ff = (ff - np.min(ff))/(np.max(ff) - np.min(ff))

    popt, pcov = curve_fit(func_fit, ww, ff.real)

    cf_fit = func_fit(ww, popt[0], popt[1], popt[2])

    ## Curve fit

    # # Gaussian weight
    # std = np.std(ww)
    # mean = np.mean(ww)
    # weight = 1/sqrt(2*pi*std**2)*np.exp(-(ww - mean)**2/2/std**2)
    # weight = weight/np.max(weight)

    # # Uniform weight
    # weight = [1 for x in ww] 

    # cf_fit = np.poly1d(np.polyfit(ww, ff.real, 2, w=weight))
    # cf_fit = cf_fit(ww)

    cf_fit_err = np.abs(1-cf_fit/ff.real)
    #cf_fit_err = ((cf_fit - ff.real) - np.min(ff.real))/(np.max(ff.real) - np.min(ff.real))

    ax_fit.plot(ww/2/pi, ff.real, color='k', linestyle='solid')
    ax_fit.plot(ww/2/pi, cf_fit, color='k', linestyle='dashed')
    ax_err.semilogy(ww/2/pi, cf_fit_err, color='k')
    #ax_err.plot(ww/2/pi, cf_fit_err, color='k')


ax_fit_ylims = ax_fit.get_ylim()
ax_fit.vlines(f_pieces, ax_fit_ylims[0], ax_fit_ylims[1], colors='k', linestyles='dashdot')
ax_fit.set_xlabel(r'Frequency [Hz]')
ax_fit.legend([r"Actual", r"$e_2 \omega^2 + e_1 \omega + e_0$"], loc='best')
ax_fit.grid(True)

ax_err_ylims = ax_err.get_ylim()
ax_err.vlines(f_pieces, ax_err_ylims[0], ax_err_ylims[1], colors='k', linestyles='dashdot')
ax_err.set_xlabel(r'Frequency [Hz]')
ax_err.set_ylabel(r'Relative error')
ax_err.legend()
ax_err.grid(True)

plt.show()
