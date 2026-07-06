"""
Plot del espectro JWST de M87 y medida de flujos del continuo con specutils
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
import scipy.constants as cst
from specutils import Spectrum, SpectralRegion

# ─────────────────────────────────────────────
# 1. Cargar espectro JWST
# ─────────────────────────────────────────────
hdu = fits.open('C:/Users/hp/Desktop/Internship IAC 2026/python-directory-almu/JWST/new/m87_full_1d_spec_corrected.fits')[1].data
wave = hdu['wave']           # micras
flux = hdu['agn_spec']*1e-3  # mJy -> Jy

# ─────────────────────────────────────────────
# 2. Medir flujos del continuo con specutils
# ─────────────────────────────────────────────
spec = Spectrum(spectral_axis=wave*u.um, flux=flux*u.Jy)

# Regiones limpias de continuo (evitando lineas de emision)
lambdas_cont = [2.8,3.3, 3.8, 4.4, 5.2, 5.7, 6.4, 7.5, 8.3, 9.8, 11.5, 13.5, 15.0, 16.2, 17.4, 20.0,21.0, 22.0,23.0, 24.0, 24.5, 25.0]  # micras
width = 0.1  # ventana de +/- 0.1 micras

print('-'*55)
print('lambda [um]    nu [Hz]         Fn [Jy]      error [Jy]')
print('-'*55)

lam_meas, nu_meas, flux_meas = [], [], []
for lam in lambdas_cont:
    mask = (wave > lam-width) & (wave < lam+width)
    fl = np.median(flux[mask])  # mediana del flujo en la ventana
    rms = np.std(flux[mask])   # error = desviacion tipica en la ventana
    nu = cst.c / (lam*1e-6)
    print(f'| {nu:>14.3e} | {fl:>14.4e} | {rms:>14.3e} | ')
    lam_meas.append(lam)
    nu_meas.append(nu)
    flux_meas.append(fl)

print('-'*55)

# ─────────────────────────────────────────────
# 3. Plot
# ─────────────────────────────────────────────
lineas = {
    'Paa':          1.875,
    '[MgIV]?':      4.487,
    'H2/[FeII]':    5.34,
    '[ArII]':       6.985,
    '[ArIII]/H2':   8.99,
    '[SIV]':        10.51,
    '[NeII]':       12.81,
    '[NeV]':        14.32,
    '[NeIII]':      15.56,
    '[FeII]':       17.94,
    '[SIII]':       18.71,
    '[OIV]':        25.89,
}

fig, ax = plt.subplots(figsize=(13, 6))

# Espectro continuo
ax.plot(wave, np.log10(flux), color='royalblue', lw=0.5, label='M87 JWST')

# Puntos de continuo medidos
ax.plot(lam_meas, np.log10(flux_meas), 'o', color='red',
        ms=6, zorder=5, label='Continuum flux measurements')

# Lineas de emision
ymin, ymax = np.nanmin(np.log10(flux)), np.nanmax(np.log10(flux))
for nombre, lam in lineas.items():
    ax.axvline(lam, ls='--', lw=0.6, color='gray', alpha=0.7)
    ax.text(lam+0.1, ymax, nombre, fontsize=7, rotation=90,
            color='gray', va='top')

ax.set_xlabel('lambda [um]', fontsize=13)
ax.set_ylabel('log(Fnu [Jy])', fontsize=13)
ax.set_xlim(2, 28)
ax.set_ylim(ymin-0.2, ymax+0.5)
ax.set_title('M87 - JWST spectrum', fontsize=14)
ax.legend(fontsize=10, loc='upper left')
plt.tight_layout()
plt.savefig('jwst_m87_spectrum.png', dpi=140)
plt.show()
