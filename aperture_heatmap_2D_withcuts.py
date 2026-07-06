import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import PowerNorm

# ============================================================
# CONFIGURACION POR IMAGEN  (lo que cambias cada vez)
# ============================================================
f = "C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/Integrated images/integrated_image_m87_nirspec_F300M_2.831-3.156um.fits"
INSTRUMENT = "NIRSpec"
BAND_LABEL = "NIRSpec F300M 2.831-3.156 um"
N_PLANOS   = 621
# CUTS de contraste en Jy (AJUSTAR por imagen para ver bien el jet)
LOW  = 1e-5
HIGH = 2.3e-4
# ============================================================
AP_PX       = 2          # radio apertura (px) - FWHM nucleo ~2px
ANN_IN_AS   = 0.4        # back interior (arcsec)
ANN_OUT_AS  = 0.7        # back exterior (arcsec)
# ============================================================

if INSTRUMENT == "NIRSpec":
    PIX_SCALE = 0.10; PIXAR_SR = 2.35040007e-13
else:
    PIX_SCALE = 0.13; PIXAR_SR = 3.97217571e-13

with fits.open(f) as h:
    names = [x.name for x in h]
    img = h["SCI"].data.copy() if "SCI" in names else h[0].data.copy()
    hdr = h["SCI"].header     if "SCI" in names else h[0].header

img = img / N_PLANOS * PIXAR_SR * 1e6   # -> Jy/pixel

yc, xc = np.unravel_index(np.nanargmax(img), img.shape)
ny, nx = img.shape
ext = [-(xc+0.5)*PIX_SCALE, (nx-xc-0.5)*PIX_SCALE,
       -(yc+0.5)*PIX_SCALE, (ny-yc-0.5)*PIX_SCALE]

ap_as = AP_PX*PIX_SCALE

fig, ax = plt.subplots(figsize=(7.5, 7.5))
# escala SQRT (PowerNorm gamma=0.5) con cuts manuales
im = ax.imshow(img, origin='lower', cmap='afmhot', extent=ext,
               norm=PowerNorm(gamma=0.5, vmin=LOW, vmax=HIGH))

# apertura (r=2px) y back (0.4-0.7")
ax.add_patch(Circle((0,0), ap_as, fill=False, ec='lime', lw=1.8))
ax.add_patch(Circle((0,0), ANN_IN_AS, fill=False, ec='cyan', ls='--', lw=1.3))
ax.add_patch(Circle((0,0), ANN_OUT_AS, fill=False, ec='cyan', ls='--', lw=1.3))
ax.plot(0, 0, '+', color='white', ms=10, mew=1.5)

# etiquetas sobre cada circulo
ax.text(0, ap_as+0.04, f'Aperture r={ap_as:.2f}"', color='lime',
        ha='center', va='bottom', fontsize=8, fontweight='bold')
ax.text(ANN_OUT_AS*0.71, ANN_OUT_AS*0.71, f'Sky {ANN_IN_AS}-{ANN_OUT_AS}"',
        color='cyan', ha='left', va='bottom', fontsize=8, fontweight='bold')

nota = (f"{BAND_LABEL}\n1 px = {PIX_SCALE}\"\n"
        f"Aperture r = {ap_as:.2f}\" ({AP_PX} px)\n"
        f"Sky = {ANN_IN_AS}-{ANN_OUT_AS}\" | sqrt scale\n"
        f"cuts: {LOW:.1e} - {HIGH:.1e} Jy")
ax.text(0.03, 0.97, nota, transform=ax.transAxes, ha='left', va='top',
        fontsize=8, family='monospace', color='white',
        bbox=dict(boxstyle='round', fc='black', ec='gray', alpha=0.6))

# limites de los ejes = extension real del cubo (todo el contenido dentro, sin recortes)
ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
ax.set_xlabel("Δx (arcsec)")
ax.set_ylabel("Δy (arcsec)")
ax.set_title(f"M87 nucleus — {BAND_LABEL}")
plt.colorbar(im, ax=ax, fraction=0.046, label="Flux per pixel (Jy)")
plt.tight_layout()
plt.show()
