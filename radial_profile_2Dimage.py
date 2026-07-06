import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# ============================================================
# CONFIGURACION POR MEDIDA  (lo unico que cambias cada vez)
# ============================================================
f = "C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/Integrated images/integrated_image_m87_miri_6.30-6.50um.fits"
INSTRUMENT = "MIRI"
N_PLANOS   = 250
BAND_LABEL = "window 6.30-6.50 um"
LAMBDA_C   = 6.40
APT_FCORR  = 2345743
APT_ERR    = 12915
# ============================================================

# parametros segun instrumento
if INSTRUMENT == "NIRSpec":
    PIX_SCALE = 0.1
    PIXAR_SR  = 2.35040007e-13
else:  # MIRI
    PIX_SCALE = 0.13
    PIXAR_SR  = 3.97217571e-13
SKY_METHOD = "median"

# apertura y anillo (en pixeles) -> se convierten a arcsec para el plot
AP_PX     = 2        # radio de apertura
ANN_IN_PX = 4        # inner sky
ANN_OUT_PX= 7        # outer sky

D_TEL = 6.5          # diametro JWST (m), para FWHM teorico

# factor de conversion a Jy (mismo para flujo y error)
def to_Jy(counts):
    return counts / N_PLANOS * PIXAR_SR * 1e6

# ---- punto fotometrico (si pasaste las medidas de APT) ----
if APT_FCORR is not None:
    F_Jy   = to_Jy(APT_FCORR)
    err_Jy = to_Jy(APT_ERR) if APT_ERR is not None else float("nan")
    print(f"PUNTO SED:  lambda = {LAMBDA_C} um")
    print(f"  F = {F_Jy:.4e} Jy  +/- {err_Jy:.2e} Jy  ({err_Jy/F_Jy*100:.1f}%)")

# ---- perfil radial ----
with fits.open(f) as h:
    img = h["SCI"].data.copy() if "SCI" in [x.name for x in h] else h[0].data.copy()

img = img / N_PLANOS * PIXAR_SR * 1e6   # -> Jy por pixel

yc, xc = np.unravel_index(np.nanargmax(img), img.shape)
ny, nx = img.shape
yy, xx = np.ogrid[:ny, :nx]
dist = np.sqrt((xx - xc)**2 + (yy - yc)**2)

sky = np.nanmedian(img[(dist > 8) & (dist <= 12)])
img_s = img - sky

# perfil en pixeles
dr = 0.3
r_edges = np.arange(0, 15 + dr, dr)
r_cent, prof = [], []
for i in range(len(r_edges) - 1):
    m = (dist >= r_edges[i]) & (dist < r_edges[i+1])
    if m.sum() > 0:
        r_cent.append((r_edges[i] + r_edges[i+1]) / 2)
        prof.append(np.nanmean(img_s[m]))
r_cent = np.array(r_cent); prof = np.array(prof)

# ajuste gaussiano (en pixeles, para el fit)
def gauss(r, A, s):
    return A * np.exp(-r**2 / (2 * s**2))
popt, _ = curve_fit(gauss, r_cent, prof, p0=[prof.max(), 1.0])
fwhm_px = 2.355 * abs(popt[1])
fwhm_arcsec = fwhm_px * PIX_SCALE

# FWHM teorico (limite de difraccion): 1.22 * lambda / D  -> arcsec
fwhm_teor = 1.025 * (LAMBDA_C * 1e-6) / D_TEL * 206265
print(f"FWHM medido  = {fwhm_px:.2f} px = {fwhm_arcsec:.3f} arcsec")
print(f"FWHM teorico = {fwhm_teor:.3f} arcsec  (1.025 lambda/D)")

# --- pasar todo a ARCSEC para el plot ---
r_cent_as = r_cent * PIX_SCALE
ap_as      = AP_PX * PIX_SCALE
ann_in_as  = ANN_IN_PX * PIX_SCALE
ann_out_as = ANN_OUT_PX * PIX_SCALE

# plot (eje X en arcsec)
fig, ax = plt.subplots(figsize=(8, 5.5))
rr_as = np.linspace(0, 15 * PIX_SCALE, 500)
rr_px = rr_as / PIX_SCALE
ax.plot(r_cent_as, prof, 'o-', ms=2.5, lw=0.6, color='black', label='Data')
ax.plot(rr_as, gauss(rr_px, *popt), 'r--', lw=1.3,
        label=f'Gaussian fit (FWHM = {fwhm_arcsec:.3f}")')
# entrada invisible: FWHM teorico en la leyenda
ax.plot([], [], ' ', label=f'FWHM theor. (1.025 λ/D) = {fwhm_teor:.3f}"')
ax.axvline(ap_as, color='magenta', lw=1.2,
           label=f'Aperture r={ap_as:.2f}"')
ax.axvspan(ann_in_as, ann_out_as, color='magenta', alpha=0.15,
           label=f'Sky annulus {ann_in_as:.2f}-{ann_out_as:.2f}"')
ax.axhline(0, color='gray', lw=0.4, ls=':')

ax.set_xlabel("Radius (arcsec)")
ax.set_ylabel("Sky-subtracted flux (Jy)")
ax.set_title(f"M87 nucleus — radial profile ({INSTRUMENT} {BAND_LABEL})")

nota = (f"1 px = {PIX_SCALE} arcsec\n"
        f"Sky = {SKY_METHOD} of annulus")
ax.text(0.985, 0.75, nota, transform=ax.transAxes, ha="right", va="top",
        fontsize=10.8, family="monospace",
        bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.8))

ax.legend(fontsize=8)
ax.grid(alpha=0.25)
plt.tight_layout()
plt.show()