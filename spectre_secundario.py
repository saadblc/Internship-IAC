import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import scipy.constants as cst

# Cargar espectro M87
hdu_m87 = fits.open('C:/Users/hp/Desktop/Internship IAC 2026/python-directory-almu/JWST/new/m87_full_1d_spec_corrected.fits')[1].data
wave_m87 = hdu_m87['wave']
flux_m87 = hdu_m87['agn_spec']*1e-3  # mJy a Jy

# Líneas de emisión (en micras)
lineas = {
    'Paα':      1.875,
    '[MgIV]?':  4.487,
    'H₂/[FeII]':5.34,
    '[ArII]':   6.985,
    '[ArIII]/H₂':8.99,
    '[SIV]':    10.51,
    '[NeII]':   12.81,
    '[NeV]':    14.32,
    '[NeIII]':  15.56,
    '[FeII]':   17.94,
    '[SIII]':   18.71,
    '[OIV]':    25.89,
}

fig, ax = plt.subplots(figsize=(12,6))

ax.plot(wave_m87, np.log10(flux_m87), color='royalblue', lw=0.5, label='M87')

# Pintar líneas verticales
for nombre, lam in lineas.items():
    ax.axvline(lam, ls='--', lw=0.6, color='gray')
    ax.text(lam+0.1, -0.6, nombre, fontsize=7, rotation=90, color='gray')

ax.set_xlabel('λ [µm]', fontsize=13)
ax.set_ylabel('log(Fν [Jy])', fontsize=13)
ax.set_xlim(2, 28)
ax.legend(fontsize=10)
plt.tight_layout()
plt.show()