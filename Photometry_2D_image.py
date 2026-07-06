import numpy as np
from astropy.io import fits

# ============================================================
# CONFIGURACION
# ============================================================
base = "C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/"
f_img = base + "integrated_image_m87_nirspec_2.00-2.10um.fits"   # tu imagen integrada (collapse)

R_NUC   = 2     # radio del circulo nuclear (px)
R_IN    = 2     # radio interior del anillo de fondo (px)
R_OUT   = 4     # radio exterior del anillo de fondo (px)

# para convertir a Jy:
N_PLANOS  = 250            # <-- EDITA: nº de planos que sumo el collapse (canal_final - canal_inicial)
PIXAR_SR  = 2.35040007e-13 # NIRSpec. (MIRI: 3.97217571e-13)
# ============================================================


# --- 1. Cargar la imagen integrada ---
with fits.open(f_img) as hdul:
    img = hdul[0].data.copy() if hdul[0].data is not None else hdul["SCI"].data.copy()
print("Imagen:", img.shape)

# --- 2. Localizar el nucleo (pixel mas brillante) ---
yc, xc = np.unravel_index(np.nanargmax(img), img.shape)
print(f"Nucleo en (y, x) = ({yc}, {xc}), valor pico = {img[yc, xc]:.1f}")

# --- 3. Mascaras circular y anular (Pitagoras) ---
ny, nx = img.shape
yy, xx = np.ogrid[:ny, :nx]
dist2 = (xx - xc)**2 + (yy - yc)**2          # distancia^2 de cada pixel al centro

m_nuc = dist2 <= R_NUC**2                      # circulo nuclear
m_ann = (dist2 > R_IN**2) & (dist2 <= R_OUT**2)  # anillo de fondo

# --- 4. Fotometria (en unidades de collapse: MJy/sr * N_planos) ---
F_nuc = np.nansum(img[m_nuc]);  N_nuc = m_nuc.sum()
F_ann = np.nansum(img[m_ann]);  N_ann = m_ann.sum()
bkg_per_pix = F_ann / N_ann                    # fondo por pixel

F_corr = F_nuc - bkg_per_pix * N_nuc           # nuclear corregido de fondo

print(f"\n--- En unidades de collapse ---")
print(f"Nucleo:  N={N_nuc} px,  F_nuc = {F_nuc:.1f}")
print(f"Anillo:  N={N_ann} px,  fondo/pix = {bkg_per_pix:.1f}")
print(f"Corregido: F_corr = {F_corr:.1f}")

# --- 5. Conversion a Jy ---
# collapse SUMO los planos -> dividir entre N_PLANOS recupera el brillo medio (MJy/sr)
# luego * PIXAR_SR (quita /sr y suma sobre los px de apertura) * 1e6 (MJy->Jy)
F_nuc_Jy  = (F_nuc  / N_PLANOS) * PIXAR_SR * 1e6
F_corr_Jy = (F_corr / N_PLANOS) * PIXAR_SR * 1e6

print(f"\n--- En Jy (N_planos={N_PLANOS}) ---")
print(f"F_nuc  = {F_nuc_Jy:.4e} Jy")
print(f"F_corr = {F_corr_Jy:.4e} Jy")