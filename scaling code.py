import numpy as np
import pandas as pd

Fn = pd.read_csv("C:\\Users\\hp\\Desktop\\Internship IAC 2026\\python-directory-almu\\high\\m87_sedhigh.dat", delimiter = '|', skiprows = 11)
Fn = Fn.apply(lambda x: pd.to_numeric(x.str.strip(), errors='coerce') if x.dtype == 'object' else x)



# Flujos
Fu = Fn[(Fn.iloc[:,2]>6.) & (Fn.iloc[:,2]<30.)]

Fu_JWST = Fu[Fu.iloc[:,4] == 1.]
Fu_other = Fu[Fu.iloc[:,4] != 1.]

Fnu_JWST = Fu_JWST.iloc[:,2]
Fnu_other = Fu_other.iloc[:,2]

JWST_Fnu = Fnu_JWST.values
other_Fnu = Fnu_other.values

mean_other = np.mean(other_Fnu)
std_other = np.std(other_Fnu)


print("Mean of other Fnu:", mean_other)
print("Standard deviation of other Fnu:", std_other)

print("JWST mean constrained:", np.mean(JWST_Fnu[(JWST_Fnu < 14.) & (JWST_Fnu > 10.5)]))
print("JWST std constrained:", np.std(JWST_Fnu[(JWST_Fnu < 14.) & (JWST_Fnu > 10.5)]))

# frecuencias
nu_JWST = Fu_JWST.iloc[:,1]
nu_other = Fu_other.iloc[:,1]

JWST_nu = nu_JWST.values
other_nu = nu_other.values
mean_other_nu = np.mean(other_nu)
std_other_nu = np.std(other_nu)

print("JWST frequencies:", JWST_nu)
print("Other frequencies:", other_nu)
print("Mean of other frequencies:", mean_other_nu)
print("Standard deviation of other frequencies:", std_other_nu)


print("SCALING FACTOR:", mean_other/np.mean(JWST_Fnu[(JWST_Fnu < 14.) & (JWST_Fnu > 10.5)]))