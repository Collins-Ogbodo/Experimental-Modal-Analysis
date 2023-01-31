from DataPreprocessing import DataPrep
from Rational_Polynomial_Fraction_Method import RFPM
from Global_Rational_Polynomial_Fraction_Method import GRFPM
from PolyMAX import PolyMAX
from PolyMaxDataPrep import PolyMaxDataPrep
from StabilizationDiagram import StabDia

#Data Preprocessing
iters = [1]
reps = [1]
test_series = "BR_AR"
frf, freq = DataPrep(iters, reps, test_series)

#Applying OMA 
sensor_name = 'EXH'
N = [i for i in range(2,30)]
nat_freqs =[]
damp_ratio =[]
order = []
frf_est = []
nat_freqs_G =[]
damp_ratio_G =[]
Order_G = []
frf_est_G = []
for i in N:
    nat_freq_, dp, order_, FRF, Freq, FRF_est_ = RFPM(frf, freq, 40, 60, sensor_name, i)
    nat_freq_, dp, order_, FRF_G, Freq, FRF_est_ = GRFPM(frf, freq, 40, 60, i)
    nat_freqs.append(nat_freq_)
    order.append(order_)
    damp_ratio.append(dp)
    frf_est.append(FRF_est_)
    nat_freqs_G.append(nat_freq_)
    Order_G.append(order_)
    damp_ratio_G.append(dp)
    frf_est_G.append(FRF_est_)


#Plot the stabilization Diagram    
plot = StabDia(nat_freqs_G, FRF_G,frf_est_G, Freq, Order_G, 'Global-RFPM')
plot = StabDia(nat_freqs, FRF,frf_est, Freq, order, sensor_name)


fRF = PolyMaxDataPrep(frf)



dr, w_n= PolyMAX(fRF, freq, 47, 55, 20)