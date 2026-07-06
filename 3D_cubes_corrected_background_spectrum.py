import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# ============================================================
# CONFIGURACION
# ============================================================
base = "C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/M87 nuclei/R02-24/"
f_nuc  = base + "nuclear_spectrum_R2_m87_nirspec_2.87-5.27micras.fits"   # espectro nuclear
f_back = base + "background_ring_R2-4_spectrum_m87_nirspec_2.87-5.27micra.fits"         # <-- EDITA: tu fits de background

# numero de pixeles de cada apertura (los metes a mano)
Nnuc  = 13     # <-- EDITA: pixeles de la apertura nuclear
Nback = 36     # <-- EDITA: pixeles de la apertura de background

# que plot quieres: 1 = los tres apilados | 0 = solo el corregido
MODO = 0

# PIXAR_SR de respaldo si el 1D no lo tiene (MIRI). Para NIRSpec: 2.35040007e-13
PIXAR_FALLBACK = 3.97217571e-13
# ============================================================


def cargar_espectro(path):
    """Carga un espectro 1D y devuelve (wave_um, flujo_Jy)."""
    with fits.open(path) as hdul:
        flujo = hdul["SCI"].data.copy()      # MJy/sr
        cab   = hdul["SCI"].header
    n = cab["NAXIS1"]
    i = np.arange(n)
    wave = cab["CRVAL1"] + (i + 1 - cab["CRPIX1"]) * cab["CDELT1"]
    pixar_sr = cab["PIXAR_SR"] if "PIXAR_SR" in cab else PIXAR_FALLBACK
    flujo_Jy = flujo * 1e6 * pixar_sr        # MJy/sr -> Jy
    return wave, flujo_Jy


# --- 1. Cargar nuclear y background ---
wave, F_nuc  = cargar_espectro(f_nuc)
_,    F_back = cargar_espectro(f_back)

# --- 2. Espectro corregido: F_nuc - (Nnuc/Nback) * F_back ---
factor = Nnuc / Nback
F_corr = F_nuc - factor * F_back
print(f"Factor Nnuc/Nback = {Nnuc}/{Nback} = {factor:.4f}")

# --- 3. Plot ---
if MODO == 1:
    # los tres apilados, uno debajo de otro
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)

    fig.suptitle("M87 nuclear - JWST/NIRspec p1293 g235 (2.87–5.27 µm)", fontsize=15, fontweight="bold")
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)   # hueco para el titulo grande
    axes[0].plot(wave, F_nuc, lw=0.6, color="darkred")
    axes[0].set_title("M87 nuclear (full)")
    axes[0].set_ylabel("Flux (Jy)")

    axes[1].plot(wave, F_back, lw=0.6, color="steelblue")
    axes[1].set_title("Background")
    axes[1].set_ylabel("Flux (Jy)")

    axes[2].plot(wave, F_corr, lw=0.6, color="black")
    axes[2].set_title("Nuclear background-corrected")
    axes[2].set_ylabel("Flux (Jy)")
    axes[2].set_xlabel(r"Wavelength ($\mu$m)")

    axes[0].text(0.98, 0.04,
        f"Aperture: circular r=2 px (R1=0, R2=2), sum\n"
        f"$Nnuc$={Nnuc}, $Nback$={Nback}, factor={factor:.3f}",
        transform=axes[0].transAxes, va="bottom", ha="right",
        fontsize=8, family="monospace",
        bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.8))
    axes[1].text(0.98, 0.04,
        f"Aperture: annular px (R1=2, R2=4), sum\n"
        f"$Nnuc$={Nnuc}, $Nback$={Nback}, factor={factor:.3f}",
        transform=axes[1].transAxes, va="bottom", ha="right",
        fontsize=8, family="monospace",
        bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.8))

    for ax in axes:
        ax.grid(alpha=0.3)

else:
    # solo el corregido
    fig, ax = plt.subplots(figsize=(11, 4.5))
    fig.suptitle("M87 nuclear - JWST/NIRspec p1293 g235 (2.87–5.27 µm)", fontsize=15, fontweight="bold")
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)   # hueco para el titulo grande
    ax.plot(wave, F_corr, lw=0.6, color="black")
    ax.set_title("M87 nuclear background-corrected")
    ax.set_xlabel(r"Wavelength ($\mu$m)")
    ax.set_ylabel("Flux (Jy)")
    ax.text(0.98, 0.04,
        f"Aperture: circular r=2 px (R1=0, R2=2), sum\n"
        f"Background subtracted (Nnuc={Nnuc}, Nback={Nback})",
        transform=ax.transAxes, va="bottom", ha="right",
        fontsize=8, family="monospace",
        bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.8))
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.show()