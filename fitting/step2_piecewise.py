from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt
import deal_vec
import piecewise_polynomial_fit as ppf
import vectfit as vf

# No. of data files
n_steps = 300

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

for n_eig in range(80, 100):

    kt_recon = np.dot(phi[:, :n_eig],  a[:n_eig, :])

    err = np.array([])
    for k in range(n_steps):
        err = np.append(err, np.sqrt(np.mean((kt_ratio[:, k]-kt_recon[:, k])**2))/(np.max(kt_ratio[:, k]) - np.min(kt_ratio[:, k])))

    ax.semilogy(np.linspace(1, n_steps, n_steps), err, label='n = {}'.format(n_eig+1))

print('Done!')

ax.set_ylabel('Normalised RMS error')
ax.set_xlabel('Step')
ax.legend()

#for k in range(len(a[:, 0])):
#f_fit, c_fit = ppf.curve_fit(a[0,:], 2, 150, 750, 300, 3)
#c_fit = [poly1d(cc) for cc in c_fit]

# ## Vector fit
# fig_fit = plt.figure(1)
# ax_fit = fig_fit.add_subplot(111)
# fig_err = plt.figure(2)
# ax_err = fig_err.add_subplot(111)
# 
# vf_fit = []
# vf_fit_err = []
# 
# ww = np.linspace(150, 750, 300)
# ff = a[1, :]
# ff = (ff - np.min(ff))/(np.max(ff) - np.min(ff))
# ax_fit.plot(ww, ff, color='k', linestyle='solid')
# for k in range(1, 6):
#     poles, residues, d, h = vf.vectfit_auto_rescale(ff, ww, n_poles=k)
#     vf_fit.append(vf.model(ww, poles, residues, d, h))
#     ax_fit.plot(ww, vf_fit[k-1].real, label="N = {}".format(k))
# 
#     vf_fit_err.append(np.abs(1-vf_fit[k-1].real/ff.real))
#     ax_err.semilogy(ww, vf_fit_err[k-1], label="N = {}".format(k))
# 
# ax_fit.legend()
# ax_err.legend()
# plt.show()

## Piecewise polynomial fit

## Global frequency range
f_0 = 150 # Hz
f_1 = 750 # Hz
w_g = np.linspace(2*pi*f_0, 2*pi*f_1, n_steps) # radians

pieces = 30
piece_len = (f_1 - f_0)/pieces

index_pieces = [0]
f_pieces = [f_0]
while True:
    if f_pieces[-1] >= f_1:
        break
    else:
        index_pieces.append(np.argmin(np.abs(w_g - (f_pieces[-1]+piece_len)*2*pi)))
        f_pieces.append(w_g[index_pieces[-1]]/2/pi)


## Figures
fig_fit = plt.figure(1)
ax_fit = fig_fit.add_subplot(111)

#fig_err = plt.figure(2)
#ax_err = fig_err.add_subplot(111)

## Fitting loop
n_f = 0
for n_ii, ii in enumerate(index_pieces[:-1]):

    if n_ii < len(index_pieces)-2:
        ii_1 = index_pieces[n_ii+1] - 1
    else:
        ii_1 = index_pieces[n_ii+1]

    ww = w_g[ii:ii_1]

    ff = a[2, ii:ii_1]
    #ff = (ff - np.min(ff))/(np.max(ff) - np.min(ff))


    ## Curve fit

    # Gaussian weight
    std = np.std(ww)
    mean = np.mean(ww)
    weight = 1/sqrt(2*pi*std**2)*np.exp(-(ww - mean)**2/2/std**2)
    weight = weight/np.max(weight)

    # Uniform weight
    weight = [1 for x in ww] 

    cf_fit = np.poly1d(np.polyfit(ww, ff.real, 2, w=weight))
    cf_fit = cf_fit(ww)

    #cf_fit_err = np.abs(1-cf_fit/ff.real)
    cf_fit_err = cf_fit - ff.real

    ax_fit.plot(ww/2/pi, ff.real, color='k', linestyle='solid')
    ax_fit.plot(ww/2/pi, cf_fit, color='k', linestyle='dashed')
    #ax_err.semilogy(ww/2/pi, cf_fit_err, color='k')
    ax_fit.plot(ww/2/pi, cf_fit_err, color='r')


ax_fit_ylims = ax_fit.get_ylim()
ax_fit.vlines(f_pieces, ax_fit_ylims[0], ax_fit_ylims[1], colors='k', linestyles='dashdot')
ax_fit.set_xlabel(r'Frequency [Hz]')
ax_fit.legend([r"Actual", r"$e_2 \omega^2 + e_1 \omega + e_0$"], loc='best')
ax_fit.grid(True)

#ax_err_ylims = ax_err.get_ylim()
#ax_err.vlines(f_pieces, ax_err_ylims[0], ax_err_ylims[1], colors='k', linestyles='dashdot')
#ax_err.set_xlabel(r'Frequency [Hz]')
#ax_err.set_ylabel(r'Relative error')
#ax_err.legend()
#ax_err.grid(True)

plt.show()
