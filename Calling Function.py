from DataPreprocessing import DataPrep
from Rational_Polynomial_Fraction_Method import RFPM
from Global_Rational_Polynomial_Fraction_Method import GRFPM
from PolyMAX import PolyMAX
from PolyMAXFFT import PolyMAXFFT
from PolyMaxData import PolyMaxDataPrep
from StabilizationDiagram import StabDia

#Data Preprocessing
iters = [1]
reps = [1]
test_series = "BR_AR"
frf, freq, coh = DataPrep(iters, reps, test_series)

#Applying OMA 
N = [i for i in range(0,15)]
min_freq = 5.0
max_freq = 55.0
#[5.0, 10.3, 13.5, 21.8, 25.52, 30.07, 39.0, 48.0, 55.0]
#%%
#RFPM parameters
sensor_name = 'EXH'
nat_freqs =[]
damp_ratio =[]
order = []
frf_est = []
#N = [1]
#OMA for multiple order  nat_freq, dam_ratio, N, FRF, Freq
for i in N:
    wn, dp, Order, FRF, Freq, FRF_est = RFPM(frf, freq, min_freq, max_freq, sensor_name, i)
    #Natural frequency
    nat_freqs.append(wn)
    #Damping Ratio
    damp_ratio.append(dp)
    #order from method
    order.append(Order)
    #Estimated FRF
    frf_est.append(FRF_est)
#Plot the stabilization Diagram    
plot = StabDia(nat_freqs, FRF,frf_est, Freq, order, sensor_name)
#%%
#GRFPM parameters
nat_freqs_G =[]
damp_ratio_G =[]
order_G = []
frf_est_G = []
#N= [16]
#OMA for multiple order  nat_freq, dam_ratio, N, FRF, Freq
for i in N:
    wn_G, dp_G, Order_G, fRF_G, _, FRF_est_G = GRFPM(frf, freq, min_freq, max_freq, i)
    #Natural frequency
    nat_freqs_G.append(wn_G)   
    #Damping Ratio
    damp_ratio_G.append(dp_G)
    #order from method
    order_G.append(Order_G)
    #Estimated FRF
    frf_est_G.append(FRF_est_G)
#Plot the stabilization Diagram    
plot = StabDia(nat_freqs_G, fRF_G,frf_est_G, Freq, order_G, 'Global-RFPM')

#%%
#Polymax
fRF = PolyMaxDataPrep(frf,['ULE_06','ULE_07'])
cOH = PolyMaxDataPrep(coh,['ULE_06','ULE_07'])
Nmin = 150
Nmax = 200
#OMA for multiple order  nat_freq, dam_ratio, N, FRF, Freq
wn_P, dp_P, Order_P, FRF_P, Freq_P = PolyMAX(fRF, freq, cOH, min_freq, max_freq, Nmin, Nmax)

#Plot the stabilization Diagram    
plot = StabDia(wn_P,FRF_P, _, Freq_P, Order_P, 'PolyMAX', 'no', 'P')

#%%

#Polymax
fRF = PolyMaxDataPrep(frf,['ULE_06','ULE_07'])
cOH = PolyMaxDataPrep(coh,['ULE_06','ULE_07'])
Nmin = 0
Nmax = 500
#OMA for multiple order  nat_freq, dam_ratio, N, FRF, Freq
wn_P, dp_P, Order_P, FRF_P, Freq_P, imax, imin = PolyMAXFFT(fRF, freq, cOH, min_freq, max_freq, Nmin, Nmax)

#Plot the stabilization Diagram    
plot = StabDia(wn_P,FRF_P, _, Freq_P, Order_P, 'PolyMAX', 'no','P')
#%%

def FreqSeg(FRF, Freq, seg, test_series,start_sensor, end_sensor):
    import matplotlib.pyplot as plt 
    import numpy as np
    import pandas as pd
    FRF = pd.DataFrame(FRF)
    # Create figure and subplot
    imin = Freq.index(0.0)
    imax = Freq.index(60)
    #converting all data to array
    Freq = np.array(Freq)
    #Frequency, Coh and FRF range
    Freq = Freq[imin:imax]
    FRF = FRF.iloc[imin:imax,start_sensor:end_sensor]
    import numpy as np
    plt.figure(figsize=(17, 8))
    plt.xlim(0,max(Freq))
    plt.grid()  
    plt.xlabel("Frequency[Hz]")
    plt.ylabel("FRF")
    for i in range(np.shape(FRF)[1]):
        FRFs = [abs(frf) for frf in FRF.iloc[:,i]]
        plt.semilogy(Freq, FRFs, label = FRF.iloc[:, i].name)
    # multiple segments
    plt.vlines(seg, ymin=0, ymax=1, colors='purple', ls='--', lw=2, label='EMA Segments')
    #for i in range(len(seg)):
       # plt.text(seg[i],seg[i+1], 'Segment'+str(i),  ha='center', va='bottom' )
    plt.title('Modal Analysis Frequency Segment-'+ test_series)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=10)
    plt.show()
    return
#%%
for i in range(0, 59, 2):
    plot = FreqSeg(frf, freq, [5.0, 10.3, 13.5, 21.8, 25.52, 30.07, 39.0, 48.0, 55.0],test_series,i, i+3)

