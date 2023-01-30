from DataPreprocessing import DataPrep
<<<<<<< Updated upstream
from Rational_Polynomial_Fraction_Method import RFPM
from StabilizationDiagram import StabDia

#Data Preprocessing
=======
from Rational Polynomial Fraction Method import Data_Prep
from Data_Preprocessing import Data_Prep
>>>>>>> Stashed changes
iters = [1]
reps = [1]
test_series = "BR_AR"
frf, freq = DataPrep(iters, reps, test_series)

#Applying OMA 
sensor_name = 'EXH'
N = [i for i in range(2,6)]
nat_freqs =[]
damp_ratio =[]
order = []
frf_est = []
for i in N:
    nat_freq_, dp, order_, FRF, Freq, FRF_est_ = RFPM(frf, freq, 30, 50, sensor_name, i)
    nat_freqs.append(nat_freq_)
    order.append(order_)
    damp_ratio.append(dp)
    frf_est.append(FRF_est_)

#Plot the stabilization Diagram    
plot = StabDia(nat_freqs, FRF,frf_est, Freq, order, sensor_name)
