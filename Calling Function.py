from DataPreprocessing import DataPrep
from Rational_Polynomial_Fraction_Method import RFPM
from Global_Rational_Polynomial_Fraction_Method import GRFPM
from PolyMAX import PolyMAX
from PolyMaxData import PolyMaxDataPrep
from StabilizationDiagram import StabDia

#Data Preprocessing
iters = [1]
reps = [1]
test_series = "BR_AR"
frf, freq = DataPrep(iters, reps, test_series)

#Applying OMA 
N = [i for i in range(2,30)]
min_freq = 20
max_freq = 60
#%%
#RFPM parameters
sensor_name = 'EXH'
nat_freqs =[]
damp_ratio =[]
order = []
frf_est = []
#OMA for multiple order  nat_freq, dam_ratio, N, FRF, Freq
for i in N:
    wn, dp, Order, FRF, Freq, FRF_est = RFPM(frf, freq, min_freq, max_freq, sensor_name, i)
    #Natural frequency
    nat_freqs.append(wn)
    #Damping Ratio
    damp_ratio.append(dp)
    #order from method
    order.append(order)
    #Estimated FRF
    frf_est.append(FRF_est)
#Plot the stabilization Diagram    
plot = StabDia(nat_freqs, FRF,frf_est, Freq, order, sensor_name, 'yes')
#%%
#GRFPM parameters
nat_freqs_G =[]
damp_ratio_G =[]
order_G = []
frf_est_G = []
#OMA for multiple order  nat_freq, dam_ratio, N, FRF, Freq
for i in N:
    wn_G, dp_G, Order_G, FRF_G, _, FRF_est_G = GRFPM(frf, freq, min_freq, max_freq, i)
    #Natural frequency
    nat_freqs_G.append(wn_G)   
    #Damping Ratio
    damp_ratio_G.append(dp_G)
    #order from method
    order_G.append(order_G)
    #Estimated FRF
    frf_est_G.append(FRF_est_G)
#Plot the stabilization Diagram    
plot = StabDia(nat_freqs_G, FRF_G,frf_est_G, Freq, Order_G, 'Global-RFPM', 'yes')
#%%
#Polymax
fRF = PolyMaxDataPrep(frf)
nat_freqs_P =[]
damp_ratio_P = []
order_P = []

#OMA for multiple order  nat_freq, dam_ratio, N, FRF, Freq
for i in N:
    wn_P, dp_P, Order_P, _, _ = PolyMAX(fRF, freq, min_freq, max_freq, i)
    #Natural frequency
    nat_freqs_P.append(wn_P)  
    #Damping Ratio
#Plot the stabilization Diagram    
plot = StabDia(nat_freqs_P, FRF_G, _, Freq, order, 'PolyMAX', 'no')


