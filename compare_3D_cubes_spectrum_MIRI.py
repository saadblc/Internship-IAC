import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# ─────────────────────────────────────────────
base = "C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/M87 nuclei/R02-24/"

f_nuc_old = base + "nuclear_spectrum_R2_m87_miri_4.90-7.65micra.fits"
f_nuc_new = base + "nuclear_spectrum_R2_m87_miri_4.90-7.65micra_new.fits"
f_bck_old = base + "background_ring_R2-4_spectrum_m87_miri_4.90-7.65micra.fits"
f_bck_new = base + "background_ring_R2-4_spectrum_m87_miri_4.90-7.65micra_new.fits"

NNUC, NBACK = 13, 36
FACTOR = NNUC / NBACK          # = 0.361
PIXAR_FALLBACK = 3.97217571e-13   # MIRI


def cargar_espectro(path):
    with fits.open(path) as hdul:
        flujo = hdul["SCI"].data.copy()
        cab   = hdul["SCI"].header
    n = cab["NAXIS1"]
    i = np.arange(n)
    wave = cab["CRVAL1"] + (i + 1 - cab["CRPIX1"]) * cab["CDELT1"]
    pixar_sr = cab["PIXAR_SR"] if "PIXAR_SR" in cab else PIXAR_FALLBACK
    return wave, flujo * 1e6 * pixar_sr        # MJy/sr -> Jy


# cargar
w_old, Fn_old = cargar_espectro(f_nuc_old)
w_new, Fn_new = cargar_espectro(f_nuc_new)
_,     Fb_old = cargar_espectro(f_bck_old)
_,     Fb_new = cargar_espectro(f_bck_new)

# corregidos
Fc_old = Fn_old - FACTOR * Fb_old
Fc_new = Fn_new - FACTOR * Fb_new

# ─────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=True)

# Panel 1: full
ax1.plot(w_old, Fn_old, lw=0.8, color="tab:red",    label="MIRI old")
ax1.plot(w_new, Fn_new, lw=0.8, color="tab:blue", label="MIRI new")
ax1.set_title("Full spectra (no background subtraction)")
ax1.set_ylabel("Flux density (Jy)")
ax1.set_ylim(0.003, 0.009)        # AJUSTA para recortar spikes
ax1.legend(fontsize=9)
ax1.grid(alpha=0.3)

# Panel 2: background-subtracted
ax2.plot(w_old, Fc_old, lw=0.8, color="tab:red",    label="MIRI old")
ax2.plot(w_new, Fc_new, lw=0.8, color="tab:blue", label="MIRI new")
ax2.set_title("Background-subtracted spectra")
ax2.set_xlabel("Wavelength (μm)")
ax2.set_ylabel("Flux density (Jy)")
ax2.set_ylim(0.002, 0.008)        # AJUSTA
ax2.legend(fontsize=9)
ax2.grid(alpha=0.3)

fig.suptitle("M87 nucleus (R2) — MIRI old vs new",
             fontsize=14, fontweight="bold")
plt.tight_layout()
plt.subplots_adjust(top=0.93)
plt.show()