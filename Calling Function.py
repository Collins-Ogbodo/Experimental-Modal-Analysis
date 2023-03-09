%reset -f
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
#%%
# ranges = ((min_freq, max_freq),(n_mode, num_ord))
# =============================================================================
# MODE-1 ((5, 13), (6, 1)),
# MODE-2 ((5, 15), (6,1)),
# MODE-3 ((15, 22), (5,1)),
# MODE-4 ((13, 21), (5,1)),
# MODE-5 ((17, 22.5), (4,1))
# MODE-6 ((18, 25), (4,1)),
# MODE-7 ((18, 25), (4,1)),
# MODE-8 ((25, 39), (5,1)) AND ((26, 39), (5,1))
# MODE-9 ((25, 39), (5,1))
# MODE-10 ((25, 39), (5,1))
# MODE-11 ((39, 43), (3,1))
# MODE-12 ((42, 47), (3,1)) and ((42, 47), (3,2))
# MODE-13 ((47, 55), (2,2))
# MODE-14 ((47, 55), (2,2))
#[5.0, 10.3, 13.5, 21.8, 25.52, 30.07, 39.0, 48.0, 55.0]
# =============================================================================
ranges = (((5, 13), (6, 1)),)

#RFPM parameters
sensor = list(frf.keys())
sensors = ['UTC_03']
for min_max_freq, n_mode_num_ord in ranges:
    print(n_mode_num_ord)
    for sensor_name in sensor:
        nat_freqs =[]
        damp_ratio =[]
        order = []
        frf_est = [] 
        #OMA for multiple order  nat_freq, dam_ratio, N, FRF, Freq
        wn, dp, Order, FRF, Freq, FRF_est = RFPM(frf, freq, min_max_freq[0], min_max_freq[1], sensor_name, n_mode_num_ord[0], n_mode_num_ord[1])
        #Natural frequency
        nat_freqs.append(wn)
        #Damping Ratio
        damp_ratio.append(dp)
        #order from method
        order.append(Order)
        #Estimated FRF
        frf_est.append(FRF_est)
        #Plot the stabilization Diagram    
        plot = StabDia(nat_freqs, FRF,frf_est, Freq, order, sensor_name, test_series, iters)
#%%
#GRFPM parameters
nat_freqs_G =[]
damp_ratio_G =[]
order_G = []
frf_est_G = []
N= [1]
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
#plot = StabDia(nat_freqs_G, fRF_G,frf_est_G, Freq, order_G, 'Global-RFPM')

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

def FreqSeg(FRF, Freq, seg, test_series,start_sensor, end_sensor, counter):
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
    plt.savefig('Results\\'+test_series+str(counter)+".png")
    plt.show()
    return
#%%
plot = FreqSeg(frf, freq, [7.086238618, 12.76555162	, 15.55088081, 17.10720227, 19.98401053, 22.74799942,	23.90120841,	27.54234413, 31.62898535, 34.9157204, 41.00430037, 43.08228857, 49.50642848, 51.85985068
],test_series,i, i+3)
#%%
counter = 0 
for i in range(0, 59, 2):
    plot = FreqSeg(frf, freq, [7.086238618, 12.76555162	, 15.55088081, 17.10720227, 19.98401053, 22.74799942,	23.90120841,	27.54234413, 31.62898535, 34.9157204, 41.00430037, 43.08228857, 49.50642848, 51.85985068
    ],test_series,i, i+3, counter)
    counter +=1

