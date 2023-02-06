import numpy as np
import json

# %% PARAMETERS

fname = 'H_simulated'
m = 1 # mass (same for all masses)
c = 0.01 # prop. damping (same for all masses)
k = 10 # stiffness (same for all masses)
n_modes = 5 # true model order
w_range = (0, 10) # frequency range (Hz)
n_w = 1000 # spectral lines
RNG_seed = 111
RMS_noise_percent = 0.0 # RMS noise level (as percentage)

# %% Helper function

def get_modal_transfer(Phi, wn2s, zwn, ws):
    return (Phi[0,:]*((-ws**2 +2j*zwn[:,None]*ws + wn2s[:, None])**-1).T)@Phi.T

# %% Generate the data

# Build parameter matricies
M = np.diag([m]*n_modes)
K = (k)*(2*M - np.diag([1]*(n_modes-1), 1) - np.diag([1]*(n_modes-1), -1))
C = c*K

# Get the modal properties
wn2s, Phi_true = np.linalg.eig(np.linalg.inv(M)@K)
wns_true = np.sqrt(wn2s)
zetas_true = np.diag(Phi_true.T@C@Phi_true) / (2*wns_true)

# Reorder modes by increasing frequency
sort_idx = np.argsort(wns_true)
wns_true = wns_true[sort_idx]
zetas_true = zetas_true[sort_idx]
Phi_true = Phi_true[:, sort_idx]
Phi_true /= np.linalg.norm(Phi_true, axis=0)*-np.sign(Phi_true[0, :])

# Get FRF
ws = np.linspace(*w_range, n_w)
H = get_modal_transfer(Phi_true, wns_true**2, zetas_true*wns_true, ws)
Ht = H

# Corrupt with noise
RNG = np.random.default_rng(seed=RNG_seed)
Hr = np.real(Ht)+(0.01*RMS_noise_percent*np.std(np.real(Ht))
                  * RNG.normal(0, 1, size=Ht.shape))
Hi = np.imag(Ht)+(0.01*RMS_noise_percent*np.std(np.imag(Ht))
                  * RNG.normal(0, 1, size=Ht.shape))
H = Hr + 1j*Hi

# save data
np.save(f'{fname}.npy',  {'FRF':H, 'ws':ws, 'wns':wns_true, 'zs':zetas_true, 'Phi':Phi_true}, allow_pickle=1)  

# load and plot the data
data = np.load(f'{fname}.npy', allow_pickle=1)[()]
import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, figsize=(10, 6))
axs[0].semilogy(data['ws'], np.abs(data['FRF']))
axs[1].plot(data['ws'], np.angle(data['FRF']))
[[ax.axvline(x, c='k', ls='--') for x in data['wns']] for ax in axs]
plt.show()