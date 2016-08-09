from __future__ import division
from __future__ import print_function

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from math import pi, sqrt

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
w_g = np.linspace(2*pi*f_0, 2*pi*f_1, 600) # radians

pieces = 2
w_g_split = np.array_split(w_g, pieces)
f_pieces = [ww[0]/2/pi for ww in w_g_split]
f_pieces.append(w_g[-1]/2/pi)
#f_pieces = [[ww[0]/2/pi,ww[-1]/2/pi] for ww in w_g_split]
#f_pieces = [item for sublist in f_pieces for item in sublist]

## Figures
fig_fit = plt.figure(1)
ax_fit = fig_fit.add_subplot(111)

fig_err = plt.figure(2)
ax_err = fig_err.add_subplot(111)

def func_fit(x, e0, e1, e2):
    return e0 + e1*x + e2/x

## Fitting loop
for w in w_g_split:

    s = 1j*w

    k02 = w**2/c**2
    kv2 = -s*rho/mu
    kh2 = -s*rho*cp/kappa

    zeta = 0.5*np.sqrt(k02/kv2) + (gamma-1)*np.sqrt(k02/kh2)

    f = zeta
    # f = s.imag*np.sqrt(s.imag)

    popt, pcov = curve_fit(func_fit, s.imag, f.real)

    cf_fit = func_fit(s.imag, popt[0], popt[1], popt[2])


    ## Curve fit

    # #std = np.std(s.imag)
    # #mean = np.mean(s.imag)
    # #weight = 1/sqrt(2*pi*std**2)*np.exp(-(s.imag - mean)**2/2/std**2)
    # #weight = weight/np.max(weight)

    # weight = [1 for x in s.imag]

    # cf_fit = np.poly1d(np.polyfit(s.imag, f.real, 2, w=weight))
    # cf_fit = cf_fit(s.imag)

    cf_fit_err = np.abs(1-cf_fit/f.real)

    ax_fit.plot(s.imag/2/pi, f.real, color='k', linestyle='solid')
    ax_fit.plot(s.imag/2/pi, cf_fit, color='k', linestyle='dashed')
    ax_err.semilogy(s.imag/2/pi, cf_fit_err, color='k')


ax_fit_ylims = ax_fit.get_ylim()
ax_fit.vlines(f_pieces, ax_fit_ylims[0], ax_fit_ylims[1], colors='k', linestyles='dashdot')
ax_fit.set_xlabel(r'Frequency [Hz]')
ax_fit.set_ylabel(r'Admittance')
ax_fit.legend([r"$\zeta$", r"$e_0 + e_1 \omega + e_2/\omega$"], loc='best')
ax_fit.grid(True)

ax_err_ylims = ax_err.get_ylim()
ax_err.vlines(f_pieces, ax_err_ylims[0], ax_err_ylims[1], colors='k', linestyles='dashdot')
ax_err.set_xlabel(r'Frequency [Hz]')
ax_err.set_ylabel(r'Relative error')
ax_err.legend()
ax_err.grid(True)

plt.show()
