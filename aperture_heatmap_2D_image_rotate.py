import numpy as np
from astropy.io import fits
from scipy.ndimage import rotate
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# ============================================================
# CONFIGURACION POR IMAGEN  (lo unico que cambias cada vez)
# ============================================================
f = "C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/Integrated images/integrated_image_m87_miri_7.20-7.40um.fits"
INSTRUMENT = "MIRI"
BAND_LABEL = "MIRI window 7.20-7.40 um"
N_PLANOS   = 250
# ============================================================
# (lo de abajo normalmente no se toca)
AP_PX, ANN_IN_PX, ANN_OUT_PX = 3, 3, 7   # apertura y anillo (px)
FOV = 2.0                      # medio campo a mostrar (arcsec)
KNOT_RMIN_PX = 5               # radio min. de busqueda del knot (excluye nucleo)
KNOT_RMAX_PX = 12              # radio max. de busqueda
# ============================================================

# --- parametros segun instrumento (se ajustan solos) ---
if INSTRUMENT == "NIRSpec":
    PIX_SCALE = 0.10
    PIXAR_SR  = 2.35040007e-13
else:  # MIRI
    PIX_SCALE = 0.13
    PIXAR_SR  = 3.97217571e-13

# --- cargar imagen integrada ---
with fits.open(f) as h:
    names = [x.name for x in h]
    img = h["SCI"].data.copy() if "SCI" in names else h[0].data.copy()
    hdr = h["SCI"].header     if "SCI" in names else h[0].header

# conversion a Jy por pixel (solo para la escala de color)
img = img / N_PLANOS * PIXAR_SR * 1e6

# --- angulo de rotacion para poner N arriba (desde la matriz CD/PC) ---
try:
    cd1, cd2 = hdr["CDELT1"], hdr["CDELT2"]
    cd = np.array([[hdr["PC1_1"]*cd1, hdr["PC1_2"]*cd1],
                   [hdr["PC2_1"]*cd2, hdr["PC2_2"]*cd2]])
    rot_deg = np.degrees(np.arctan2(-cd[0,1], cd[1,1]))
except KeyError:
    rot_deg = 0.0

# rotar a N-up (la FOTOMETRIA se hace en la imagen ORIGINAL, no en esta)
img = rotate(np.nan_to_num(img), angle=rot_deg, reshape=True, order=1, cval=0)

# --- centro (nucleo) y distancias ---
yc, xc = np.unravel_index(np.nanargmax(img), img.shape)
ny, nx = img.shape
yy, xx = np.ogrid[:ny, :nx]
dist_px = np.sqrt((xx-xc)**2 + (yy-yc)**2)

# --- localizar el knot del jet (pico fuera del nucleo) ---
mask = (dist_px > KNOT_RMIN_PX) & (dist_px < KNOT_RMAX_PX)
yk, xk = np.unravel_index(np.argmax(np.where(mask, img, -np.inf)), img.shape)
d_knot_as   = np.sqrt((xk-xc)**2 + (yk-yc)**2) * PIX_SCALE
kx_as, ky_as = (xk-xc)*PIX_SCALE, (yk-yc)*PIX_SCALE
print(f"knot a d = {d_knot_as:.3f} arcsec del nucleo")

# --- radios en arcsec ---
ap_as = AP_PX*PIX_SCALE
ai    = ANN_IN_PX*PIX_SCALE
ao    = ANN_OUT_PX*PIX_SCALE

# --- extent y contraste ---
ext = [-(xc+0.5)*PIX_SCALE, (nx-xc-0.5)*PIX_SCALE,
       -(yc+0.5)*PIX_SCALE, (ny-yc-0.5)*PIX_SCALE]
vmin, vmax = np.nanpercentile(img[img > 0], [40, 99.5])

# --- plot ---
fig, ax = plt.subplots(figsize=(7.5, 7.5))
im = ax.imshow(img, origin='lower', cmap='afmhot', extent=ext, vmin=vmin, vmax=vmax)

# circulos: apertura (verde) y anillo (cian)
ax.add_patch(Circle((0,0), ap_as, fill=False, ec='lime', lw=1.8))
ax.add_patch(Circle((0,0), ai,    fill=False, ec='cyan', ls='--', lw=1.3))
ax.add_patch(Circle((0,0), ao,    fill=False, ec='cyan', ls='--', lw=1.3))
ax.plot(0, 0, '+', color='white', ms=10, mew=1.5)

# etiquetas directamente sobre cada circulo
ax.text(0, ap_as+0.05, f'Aperture r={ap_as:.2f}"', color='lime',
        ha='center', va='bottom', fontsize=8, fontweight='bold')
ax.text(ao*0.71, ao*0.71, f'Sky annulus\n{ai:.2f}-{ao:.2f}"', color='cyan',
        ha='left', va='bottom', fontsize=8, fontweight='bold')

# knot del jet + distancia + linea
ax.plot(kx_as, ky_as, 'o', mfc='none', mec='lime', ms=20, mew=2)
ax.plot([0, kx_as], [0, ky_as], '--', color='lime', lw=0.8)
ax.annotate(f'jet knot d={d_knot_as:.2f}"', (kx_as, ky_as), color='lime',
            fontsize=8, xytext=(0,-26), textcoords='offset points', ha='center')

# barra de escala 0.5"
x0, y0 = -FOV+0.3, -FOV+0.3
ax.plot([x0, x0+0.5], [y0, y0], 'w-', lw=2.5)
ax.text(x0+0.25, y0+0.08, '0.5"', color='white', ha='center', fontsize=9)

# flechas N-E (imagen ya rotada: N arriba, E izquierda)
ox, oy, L = FOV-0.4, -FOV+0.4, 0.4
ax.annotate('', xy=(ox, oy+L), xytext=(ox, oy), arrowprops=dict(color='yellow', width=1.2, headwidth=6))
ax.text(ox, oy+1.4*L, 'N', color='yellow', ha='center', fontsize=10)
ax.annotate('', xy=(ox-L, oy), xytext=(ox, oy), arrowprops=dict(color='yellow', width=1.2, headwidth=6))
ax.text(ox-1.4*L, oy, 'E', color='yellow', va='center', fontsize=10)

# caja de anotacion
nota = (f"{BAND_LABEL}\n"
        f"1 px = {PIX_SCALE}\"\n"
        f"Aperture r = {ap_as:.2f}\" ({AP_PX} px)\n"
        f"Sky annulus = {ai:.2f}-{ao:.2f}\"\n"
        f"knot d = {d_knot_as:.2f}\"")
ax.text(0.03, 0.97, nota, transform=ax.transAxes, ha='left', va='top',
        fontsize=8, family='monospace', color='white',
        bbox=dict(boxstyle='round', fc='black', ec='gray', alpha=0.6))

ax.set_xlim(-FOV, FOV); ax.set_ylim(-FOV, FOV)
ax.set_xlabel("Δ R.A. (arcsec)")
ax.set_ylabel("Δ Dec. (arcsec)")
ax.set_title(f"M87 nucleus — {BAND_LABEL} (N-up, E-left)")
plt.colorbar(im, ax=ax, fraction=0.046, label="Flux per pixel (Jy)")
plt.tight_layout()
plt.show()