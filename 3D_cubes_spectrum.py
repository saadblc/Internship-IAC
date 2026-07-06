import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# --- 1. Cargar el espectro ---
f = "C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/M87 nuclei/nuclear_spectrum_m87_miri_4.90-7.65micra.fits"
with fits.open(f) as hdul:
    flujo = hdul["SCI"].data.copy()   # array 1D de intensidades (MJy/sr)
    cab   = hdul["SCI"].header

# --- 2. Ver la cabecera completa del SCI ---
print(repr(cab))   # vuelca todas las claves; busca PIXAR_SR aqui

# --- 3. Reconstruir el eje de longitud de onda ---
n = cab["NAXIS1"]
i = np.arange(n)
wave = cab["CRVAL1"] + (i + 1 - cab["CRPIX1"]) * cab["CDELT1"]

# --- 4. Conversion MJy/sr -> Jy ---
if "PIXAR_SR" in cab:
    pixar_sr = cab["PIXAR_SR"]
    print(f"\nPIXAR_SR encontrado en el archivo 1D: {pixar_sr}")
else:
    # si el 1D no lo guardo, lo lees del CUBO original de NIRSpec G395
    cubo_orig = "C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/M87 nuclei/m87_nirspec_1150_p1293_g395.fits"
    with fits.open(cubo_orig) as hc:
        pixar_sr = hc["SCI"].header["PIXAR_SR"]
    print(f"\nPIXAR_SR NO estaba en el 1D; leido del cubo original: {pixar_sr}")

flujo_Jy = flujo * 1e6 * pixar_sr   # MJy/sr -> Jy (suma de la apertura)

# --- 5. Dibujar en Jy ---
fig, ax = plt.subplots(figsize=(11, 4.5))
ax.plot(wave, flujo_Jy, lw=0.6, color="darkred")
ax.set_xlabel(r"Wavelength ($\mu$m)")
ax.set_ylabel("Flux (Jy)")
ax.text(0.98, 0.04, "Aperture: circular r=2 px (R1=0, R2=2), sum\nBackground: no subtraction",
transform=ax.transAxes, va="bottom", ha="right", fontsize=8, family="monospace",
bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.8))
ax.set_title("MIRI background")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()