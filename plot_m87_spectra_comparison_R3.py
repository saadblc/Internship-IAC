"""
Comparacion de espectros nucleares M87 (set R3 / R4-7):
   (nuclear R3  -  background R4-7)  convertido a mJy   vs   aperture_spec

Todos los parametros (CRVAL3, CDELT3, PIXAR_SR) se LEEN del header de cada
FITS. Cada canal MIRI tiene su propio PIXAR_SR (ch1 3.97e-13, ch2 6.79e-13, ch3/ch4 9.40e-13) y el codigo
lo toma de cada archivo automaticamente.

Procedimiento por tramo:
  1. wave = CRVAL3 + CDELT3 * arange(n)      (reconstruido del header)
  2. sub  = (nuclear - background)           (resta canal a canal)
  3. sub_mJy = sub * PIXAR_SR * 1e9          (MJy/sr -> Jy -> mJy)
  4. se superpone con 'aperture_spec' del espectro 1D completo.

"""
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# ============================================================
# RUTAS (ordenador de trabajo)
# ============================================================
UP        = r"C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/M87 nuclei/R03-47/"
FULL_SPEC = r"C:/Users/hp/Desktop/Internship IAC 2026/python-directory-almu/JWST/new/m87_full_1d_spec_corrected.fits"

# (label, nuclear R3, background R4-7, color)
PAIRS = [
    ("NIRSpec g235", "nuclear_spectrum_R3_m87_nirspec_g235.fits",
                     "background_spectrum_R4-7_m87_nirspec_g235.fits", "tab:blue"),
    ("NIRSpec g395", "nuclear_spectrum_R3_m87_nirspec_g395.fits",
                     "background_spectrum_R4-7_m87_nirspec_g395.fits", "tab:orange"),
    ("MIRI ch1",     "nuclear_spectrum_R3_m87_miri_ch1.fits",
                     "background_spectrum_R4-7_m87_miri_ch1.fits", "tab:green"),
    ("MIRI ch2",     "nuclear_spectrum_R3_m87_miri_ch2.fits",
                     "background_spectrum_R4-7_m87_miri_ch2.fits", "tab:red"),
    ("MIRI ch3",     "nuclear_spectrum_R3_m87_miri_ch3.fits",
                     "background_spectrum_R4-7_m87_miri_ch3.fits", "tab:purple"),
    ("MIRI ch4",     "nuclear_spectrum_R3_m87_miri_ch4.fits",
                     "background_spectrum_R4-7_m87_miri_ch4.fits", "tab:brown"),
]

# ============================================================
def getkey(hdr, *names):
    """Devuelve el primer keyword que exista (p.ej. CRVAL3 o CRVAL1)."""
    for n in names:
        if n in hdr:
            return hdr[n]
    raise KeyError(f"Ninguno de {names} esta en el header")

def load_1d(path):
    """Carga espectro 1D + header del HDU primario."""
    with fits.open(path) as h:
        return h[0].data.astype(float), h[0].header

# ------------------------------------------------------------
# espectro completo de referencia
with fits.open(FULL_SPEC) as h:
    full_wave = h[1].data["wave"]
    full_ap   = h[1].data["aperture_spec"]      # ya en mJy

# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(14, 6))

for label, fn_nuc, fn_bkg, color in PAIRS:
    nuc, hn = load_1d(UP + fn_nuc)
    bkg, hb = load_1d(UP + fn_bkg)

    crval = getkey(hn, "CRVAL3", "CRVAL1")
    cdelt = getkey(hn, "CDELT3", "CDELT1")
    pixar = getkey(hn, "PIXAR_SR")

    wave = crval + cdelt * np.arange(len(nuc))
    sub_mjy = (nuc) * pixar * 1e9          # MJy/sr -> mJy

    ax.plot(wave, sub_mjy, lw=0.5, color=color,
            label=f"{label} (R3)")

ax.plot(full_wave, full_ap, lw=1.0, color="k", alpha=0.7,
        label="aperture_spec (full 1D)")

ax.set_xlabel("wavelength (um)")
ax.set_ylabel("flux (mJy)")
ax.set_title("M87 nuclear: (R3 nuclear - R4-7 background) vs aperture_spec")
ax.set_xlim(1.5, 28)
ax.set_ylim(-10,50)
ax.grid(alpha=0.2)
ax.legend(fontsize=7, ncol=2)
plt.tight_layout()
plt.savefig("m87_spectra_comparison_R3.png", dpi=130, facecolor="white")
plt.show()
