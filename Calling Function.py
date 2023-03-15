%reset -f
#%%
from DataPreprocessing import DataPrep
from Rational_Polynomial_Fraction_Method import RFPM
from Global_Rational_Polynomial_Fraction_Method import GRFPM
from PolyMAX import PolyMAX
from PolyMAXFFT import PolyMAXFFT
from PolyMaxData import PolyMaxDataPrep
from StabilizationDiagram import StabDia
from collections import defaultdict
#Data Preprocessing
Iter = [1, 2]
reps = [1]
test_series = "BR_AR"
output = {}
BR_AMP_levels = [0.4, 0.8, 1.2, 1.6, 2]
DS_AMP_levels = [0.4, 1.2, 2, 0.4, 1.2, 2, 0.4, 1.2, 2]
BR_DMG_levels = [0] * 5
DS_DMG_levels = [1, 1, 1, 2, 2, 2, 3, 3, 3]
for iters, AMP_levels, DMG_levels in zip(Iter, BR_AMP_levels, BR_DMG_levels): 
    
    frf, freq = DataPrep([iters], reps, test_series)
    ranges = (
    ((5, 9), 1),
    ((12, 14), 1),
    ((15, 19), 1),
    ((22, 24), 2),
    ((26, 30), 1),
    ((30.5, 31), 1),
    ((35, 37), 1),
    ((40.5, 44), 2),
    ((48.5, 54), 2),
    ((86, 90), 1),
    ((92, 100), 1),
    ((112, 118.5), 1),
    ((119.5, 122), 1),
    ((122, 125), 1),
    ((135, 138), 1),
    ((154, 162), 1),)
    num_ord = 6
    import numpy as np
    #Experiment Description

    output[test_series+"_"+str(iters)] = {"AMP_level": AMP_levels, "DMG_level": DMG_levels, "Data":[]}
    #RFPM parameters
    for (min_freq, max_freq), n_mode in ranges:
        wns = []
        wns1 = []
        for sensor_name in list(frf.keys()):
            nat_freqs =[]
            damp_ratio =[] 
            #RFPM
            wn,_,_  = RFPM(frf, freq, min_freq, max_freq, sensor_name, n_mode, num_ord)
            #Natural frequency
            if n_mode ==2:
                rounded_numbers = [num//1 for num in wn]
                num_dict = defaultdict(list)
                for num in wn:
                    num_dict[num//1].append(num)
                modes = sorted(num_dict.values(), key=len, reverse=True)[:n_mode]
                wns.append(np.mean(modes[0]))
                wns1.append(np.mean(modes[1]))
            else:
                wns.append(np.mean(wn))
        if n_mode ==2:
            output[test_series+"_"+str(iters)]["Data"].append(np.mean(wns))
            output[test_series+"_"+str(iters)]["Data"].append(np.mean(wns1))
        else:
            output[test_series+"_"+str(iters)]["Data"].append(np.mean(wns))
            

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
plot = FreqSeg(frf, freq, [7.086238618, 12.76555162	, 15.55088081, 17.10720227, 19.98401053, 22.74799942,	23.90120841,	27.54234413, 31.62898535, 34.9157204, 41.00430037, 43.08228857, 49.50642848, 51.85985068],test_series,i, i+3)
#%%
counter = 0 
for i in range(0, 59, 2):
    plot = FreqSeg(frf, freq, [5, 9, 12, 14, 22, 24, 26, 30, 30.5, 31, 35,37,40.5, 44, 48.5, 54],test_series,i, i+3, counter)
    counter +=1

