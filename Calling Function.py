import DataPreprocessing



from Data_Preprocessing import Data_Prep
from Data_Preprocessing import Data_Prep
iters = [1]
reps = [1]
test_series = "BR_AR"

FRF, Freq = Data_Prep(iters, reps, test_series)



sensor = 'EXH'
N = [i for i in range(2,10)]
nat_freqs =[]
damp_ratio =[]
order = []
frf_est = []
for i in N:
    nat_freq_, dp, order_, FRF, Freq, FRF_est_ = RFPM_alg(sensor_frf_mean, sensor_frf_freq_mean, 20, 50, sensor, i)
    nat_freqs.append(nat_freq_)
    order.append(order_)
    damp_ratio.append(dp)
    frf_est.append(FRF_est_)
    
plot = Stabilization_Dia(FRF,frf_est, Freq, order, sensor)