from astropy.io import fits

f = "C:/Users/hp/Desktop/Internship IAC 2026/IMAGES/M87 nuclei/3Dcubos originales/m87_miri_1151_p1293_ch1.fits"
with fits.open(f) as h:
    hdr = h["SCI"].header          # OJO: header del SCI, no del primary
    print("CDELT1 =", hdr["CDELT1"], "deg")
    print("pixel  =", abs(hdr["CDELT1"])*3600, "arcsec")
    print("PIXAR_A2 =", hdr["PIXAR_A2"], "arcsec^2")
    print("lado =", hdr["PIXAR_A2"]**0.5, "arcsec")