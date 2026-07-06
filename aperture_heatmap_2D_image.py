import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# ============================================================
# CONFIGURACION POR IMAGEN (lo que cambias cada vez)
# ============================================================
f = "C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/Integrated images/integrated_image_m87_miri_7.20-7.40um.fits"
INSTRUMENT = "MIRI"
BAND_LABEL = "MIRI window 7.20-7.40 um"
N_PLANOS   = 250

AP_PX, ANN_IN_PX, ANN_OUT_PX = 3, 3, 7   # apertura y anillo (px)
JET_DIST = 1.2             # arcsec, distancia del knot al nucleo
JET_PA   = 230                # angulo posicion del jet (grados, N->E)
FOV      = 2.0                # medio campo a mostrar (arcsec)
# ============================================================

PIX_SCALE = 0.1 if INSTRUMENT == "NIRSpec" else 0.13
PIXAR_SR  = 2.35040007e-13 if INSTRUMENT == "NIRSpec" else 3.97217571e-13

with fits.open(f) as h:
    img = h["SCI"].data.copy() if "SCI" in [x.name for x in h] else h[0].data.copy()

# conversion a Jy por pixel (solo para la escala de color)
img = img / N_PLANOS * PIXAR_SR * 1e6

# centro y extent en arcsec (0,0 = nucleo)
yc, xc = np.unravel_index(np.nanargmax(img), img.shape)
ny, nx = img.shape
ext = [-(xc+0.5)*PIX_SCALE, (nx-xc-0.5)*PIX_SCALE,
       -(yc+0.5)*PIX_SCALE, (ny-yc-0.5)*PIX_SCALE]

vmin, vmax = np.nanpercentile(img, [40, 99.5])

fig, ax = plt.subplots(figsize=(7.5, 7.5))
im = ax.imshow(img, origin='lower', cmap='afmhot', extent=ext,
               vmin=vmin, vmax=vmax)

# aperturas (se dibujan a partir de AP_PX/ANN_*, se actualizan solas)
ax.add_patch(Circle((0,0), AP_PX*PIX_SCALE, fill=False, ec='lime', lw=1.8))
ax.add_patch(Circle((0,0), ANN_IN_PX*PIX_SCALE, fill=False, ec='cyan', ls='--', lw=1.3))
ax.add_patch(Circle((0,0), ANN_OUT_PX*PIX_SCALE, fill=False, ec='cyan', ls='--', lw=1.3))
ax.plot(0, 0, '+', color='white', ms=10, mew=1.5)

# marcador del knot del jet (texto debajo)
ang = np.radians(JET_PA)
jx, jy = JET_DIST*np.sin(ang), JET_DIST*np.cos(ang)
ax.plot(jx, jy, 'o', mfc='none', mec='lime', ms=26, mew=2.0)
ax.annotate("jet knot (HST-1)", (jx, jy), color='lime', fontsize=8,
            xytext=(0,-25), textcoords='offset points', ha='center')

# barra de escala 0.5"
x0, y0 = -FOV+0.3, -FOV+0.3
ax.plot([x0, x0+0.5], [y0, y0], 'w-', lw=2.5)
ax.text(x0+0.25, y0+0.08, '0.5"', color='white', ha='center', fontsize=9)

# flecha N-E
ax.annotate('', xy=(FOV-0.3, -FOV+0.7), xytext=(FOV-0.3, -FOV+0.3),
            arrowprops=dict(color='yellow', width=1.2, headwidth=6))
ax.text(FOV-0.3, -FOV+0.8, 'N', color='yellow', ha='center', fontsize=9)
ax.annotate('', xy=(FOV-0.7, -FOV+0.3), xytext=(FOV-0.3, -FOV+0.3),
            arrowprops=dict(color='yellow', width=1.2, headwidth=6))
ax.text(FOV-0.8, -FOV+0.3, 'E', color='yellow', va='center', fontsize=9)

# caja de anotacion (apertura/anillo se actualizan solos)
nota = (f"{BAND_LABEL}\n"
        f"1 px = {PIX_SCALE}\"\n"
        f"Aperture r = {AP_PX*PIX_SCALE:.2f}\" ({AP_PX} px)\n"
        f"Sky annulus = {ANN_IN_PX*PIX_SCALE:.2f}-{ANN_OUT_PX*PIX_SCALE:.2f}\"")
ax.text(0.03, 0.97, nota, transform=ax.transAxes, ha='left', va='top',
        fontsize=8, family='monospace', color='white',
        bbox=dict(boxstyle='round', fc='black', ec='gray', alpha=0.6))

ax.set_xlim(-FOV, FOV); ax.set_ylim(-FOV, FOV)
ax.set_xlabel("Δ R.A. (arcsec)")
ax.set_ylabel("Δ Dec. (arcsec)")
ax.set_title(f"M87 nucleus — {BAND_LABEL}")
plt.colorbar(im, ax=ax, fraction=0.046, label="Flux per pixel (Jy)")
plt.tight_layout()
plt.show()