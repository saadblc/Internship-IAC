"""
Comparación de agn_spec vs aperture_spec del espectro JWST de M87
y medida de flujos del continuo
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.constants as cst

# ─────────────────────────────────────────────
# 1. Cargar espectro JWST - ambas columnas
# ─────────────────────────────────────────────
hdu = fits.open('C:/Users/hp/Desktop/Internship IAC 2026/python-directory-almu/JWST/new/m87_full_1d_spec_corrected.fits')[1].data
wave          = hdu['wave']                   # micras
flux_agn      = hdu['agn_spec']        # mJy -> Jy
flux_aperture = hdu['aperture_spec']   # mJy -> Jy

mask_all = np.isfinite(flux_agn) & np.isfinite(flux_aperture)

# ─────────────────────────────────────────────
# 2. Medir flujos del continuo con ventanas
# ─────────────────────────────────────────────
lambdas_cont = [1.6, 2.1, 2.6, 2.8, 3.3, 3.8, 4.4, 5.2, 5.7, 6.4, 7.5, 8.3, 9.8,
                11.5, 13.5, 15.0, 16.2, 17.4, 20.0, 21.0, 22.0,
                23.0, 24.0, 24.5, 25.0]
width = 0.1

print('-'*75)
print(f'{"nu [Hz]":>14}  {"agn-spec [Jy]":>14}  {"error [Jy]":>12}')
print('-'*75)

lam_meas, nu_meas, flux_meas_agn, flux_meas_ap, err_meas_agn, err_meas_ap = [], [], [], [], [], []

for lam in lambdas_cont:
    mask = (wave > lam - width) & (wave < lam + width)
    if np.sum(mask) == 0:
        continue
    fl_agn = np.nanmedian(flux_agn[mask])
    fl_ap  = np.nanmedian(flux_aperture[mask])
    rms_agn   = np.nanstd(flux_agn[mask])
    rms_ap    = np.nanstd(flux_aperture[mask])
    nu     = cst.c / (lam * 1e-6)
    print(f'| {nu:>14.3e} | {fl_ap:>14.4e} |  {rms_ap:>12.3e} |')
    lam_meas.append(lam)
    nu_meas.append(nu)
    flux_meas_agn.append(fl_agn)
    flux_meas_ap.append(fl_ap)
    err_meas_agn.append(rms_agn)
    err_meas_ap.append(rms_ap)

print('-'*75)

lam_meas      = np.array(lam_meas)
flux_meas_agn = np.array(flux_meas_agn)
flux_meas_ap  = np.array(flux_meas_ap)
err_meas_agn  = np.array(err_meas_agn)
err_meas_ap   = np.array(err_meas_ap)

# ─────────────────────────────────────────────
# 3. Líneas de emisión
# ─────────────────────────────────────────────
lineas = {
    'Paα':        1.875,
    '[MgIV]?':    4.487,
    'H2/[FeII]':  5.34,
    '[ArII]':     6.985,
    '[ArIII]/H2': 8.99,
    '[SIV]':      10.51,
    '[NeII]':     12.81,
    '[NeV]':      14.32,
    '[NeIII]':    15.56,
    '[FeII]':     17.94,
    '[SIII]':     18.71,
    '[OIV]':      25.89,
}

# ─────────────────────────────────────────────
# 4. Plot único
# ─────────────────────────────────────────────
plt.rcParams['figure.dpi'] = 140
fig, ax = plt.subplots(figsize=(13, 6))

# Espectros continuos
ax.plot(wave[mask_all], np.log10(flux_agn[mask_all]),
        color='royalblue', lw=0.5, label='agn_spec', alpha=0.8)
ax.plot(wave[mask_all], np.log10(flux_aperture[mask_all]),
        color='orange', lw=0.5, label='aperture_spec', alpha=0.8)

# Puntos medidos
ax.plot(lam_meas, np.log10(flux_meas_agn), 'o',
        color='blue', ms=6, zorder=5, label='agn puntos medidos')
ax.plot(lam_meas, np.log10(flux_meas_ap), 's',
        color='darkorange', ms=6, zorder=5, label='aperture puntos medidos')

# Líneas de emisión
ymin = np.nanmin(np.log10(flux_agn[mask_all]))
ymax = np.nanmax(np.log10(flux_agn[mask_all]))
for nombre, lam in lineas.items():
    ax.axvline(lam, ls='--', lw=0.6, color='gray', alpha=0.7)
    ax.text(lam + 0.1, ymax, nombre, fontsize=7, rotation=90,
            color='gray', va='top')

ax.set_xlabel('λ [μm]', fontsize=12)
ax.set_ylabel('log(Fν [Jy])', fontsize=12)
ax.set_xlim(1.5, 28)
ax.set_ylim(ymin - 0.2, ymax + 0.5)
ax.set_title('M87 - JWST: agn_spec vs aperture_spec', fontsize=13)
ax.legend(fontsize=9, loc='upper left')

plt.tight_layout()
plt.show()

