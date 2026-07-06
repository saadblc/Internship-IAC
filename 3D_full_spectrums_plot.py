import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# ─────────────────────────────────────────────
base = "C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/M87 nuclei/R02-24/"

# nucleares (circulo r=2)
f_nuc = {
    "nir1": base + "nuclear_spectrum_R2_m87_nirspec_1.66-3.17micras.fits",
    "nir2": base + "nuclear_spectrum_R2_m87_nirspec_2.87-5.27micras.fits",
    "miri": base + "nuclear_spectrum_R2_m87_miri_4.90-7.65micra.fits",
}
# fondo (anillo r=2-4)
f_bck = {
    "nir1": base + "background_ring_R2-4_spectrum_m87_nirspec_1.66-3.17micra.fits",
    "nir2": base + "background_ring_R2-4_spectrum_m87_nirspec_2.87-5.27micra.fits",
    "miri": base + "background_ring_R2-4_spectrum_m87_miri_4.90-7.65micra.fits",
}

# numero de pixeles de cada apertura
NNUC  = 13     # circulo r=2
NBACK = 36     # anillo r=2-4
FACTOR = NNUC / NBACK   # = 0.361

PIXAR_FALLBACK = 3.97217571e-13   # MIRI (NIRSpec: 2.35040007e-13)


def cargar_espectro(path):
    """Carga un espectro 1D y devuelve (wave_um, flujo_Jy)."""
    with fits.open(path) as hdul:
        flujo = hdul["SCI"].data.copy()
        cab   = hdul["SCI"].header
    n = cab["NAXIS1"]
    i = np.arange(n)
    wave = cab["CRVAL1"] + (i + 1 - cab["CRPIX1"]) * cab["CDELT1"]
    pixar_sr = cab["PIXAR_SR"] if "PIXAR_SR" in cab else PIXAR_FALLBACK
    return wave, flujo * 1e6 * pixar_sr        # MJy/sr -> Jy


# cargar nucleares
wave_nir1, F_nir1 = cargar_espectro(f_nuc["nir1"])
wave_nir2, F_nir2 = cargar_espectro(f_nuc["nir2"])
wave_miri, F_miri = cargar_espectro(f_nuc["miri"])

# cargar fondos
_, B_nir1 = cargar_espectro(f_bck["nir1"])
_, B_nir2 = cargar_espectro(f_bck["nir2"])
_, B_miri = cargar_espectro(f_bck["miri"])

# espectros corregidos: F_nuc - FACTOR * F_back
C_nir1 = F_nir1 - FACTOR * B_nir1
C_nir2 = F_nir2 - FACTOR * B_nir2
C_miri = F_miri - FACTOR * B_miri

# ─────────────────────────────────────────────
plt.figure(figsize=(12, 6))

# originales (de colores)
plt.plot(wave_nir1, F_nir1, lw=0.8, label="NIRSpec 1.66–3.17 μm")
plt.plot(wave_nir2, F_nir2, lw=0.8, label="NIRSpec 2.87–5.27 μm")
plt.plot(wave_miri, F_miri, lw=0.8, label="MIRI 4.90–7.65 μm")

# corregidos (todos en negro)
plt.plot(wave_nir1, C_nir1, lw=0.8, color="black")
plt.plot(wave_nir2, C_nir2, lw=0.8, color="black")
plt.plot(wave_miri, C_miri, lw=0.8, color="black",
         label="Background-subtracted")

plt.xlabel("Wavelength (μm)")
plt.ylabel("Flux density (Jy)")
plt.title("M87 nucleus (R2) — JWST spectra, raw vs background-subtracted")
plt.ylim(0, 0.006)          # AJUSTA para recortar spikes
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
