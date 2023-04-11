from DataPreprocessing import DataPrep
from Rational_Polynomial_Fraction_Method import RFPM
from Global_Rational_Polynomial_Fraction_Method import GRFPM
from PolyMAX import PolyMAX
from PolyMAXFFT import PolyMAXFFT
from PolyMaxData import PolyMaxDataPrep
from StabilizationDiagram import StabDia

#Data Preprocessing
iters = [1]
reps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
test_series = "BR_AR"
frf, freq, coh = DataPrep(iters, reps, test_series)

#Applying OMA 
min_freq = 5
max_freq = 20
#%%
#Polymax
fRF = PolyMaxDataPrep(frf,['ULE_06','ULE_07'])
cOH = PolyMaxDataPrep(coh,['ULE_06','ULE_07'])
Nmin = 20
Nmax = 80
#OMA for multiple order  nat_freq, dam_ratio, N, FRF, Freq
wn_P, dp_P, Order_P, FRF_P, Freq_P = PolyMAX(fRF, freq, cOH, min_freq, max_freq, Nmin, Nmax)

#Plot the stabilization Diagram    
plot = StabDia(frf, wn_P,FRF_P, _ , Freq_P, Order_P, 'PolyMAX', 'no', 'P')

#%%

#Polymax
fRF = PolyMaxDataPrep(frf,['ULE_06','ULE_07'])
cOH = PolyMaxDataPrep(coh,['ULE_06','ULE_07'])
Nmin = 0
Nmax = 500
#OMA for multiple order  nat_freq, dam_ratio, N, FRF, Freq
wn_P, dp_P, Order_P, FRF_P, Freq_P, imax, imin = PolyMAXFFT(fRF, freq, cOH, min_freq, max_freq, Nmin, Nmax)

#Plot the stabilization Diagram    
plot = StabDia(frf, wn_P,FRF_P, _ , Freq_P, Order_P, 'PolyMAX', 'no','P')