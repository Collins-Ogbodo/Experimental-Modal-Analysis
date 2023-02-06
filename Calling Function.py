from DataPreprocessing import DataPrep
from Rational_Polynomial_Fraction_Method import RFPM
from Global_Rational_Polynomial_Fraction_Method import GRFPM
from PolyMAX import PolyMAX
from PolyMaxData import PolyMaxDataPrep
from StabilizationDiagram import StabDia

#Data Preprocessing
#iters = [1]
#reps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
#test_series = "BR_AR"
#frf, freq, coh = DataPrep(iters, reps, test_series)
#%%
frf = {}
for i in range(5):
    frf[i] = H[:,i]
 #%%   
freq = [i for i in ws]
#Applying OMA 
N = [i for i in range(9,20)]
min_freq = 0.0
max_freq = 7.007007007007007
#[5.0, 10.3, 13.5, 21.8, 25.52, 30.07, 39.0, 48.0, 55.0]
#%%
#RFPM parameters
sensor_name =2
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
fRF = PolyMaxDataPrep(frf)
cOH = 1 #PolyMaxDataPrep(coh)
nat_freqs_P =[]
damp_ratio_P = []
order_P = []

#OMA for multiple order  nat_freq, dam_ratio, N, FRF, Freq
for i in N:
    wn_P, dp_P, Order_P, _, _ = PolyMAX(fRF,cOH, freq, min_freq, max_freq, i)
    #Natural frequency
    nat_freqs_P.append(wn_P)  
    #Damping Ratio
    damp_ratio_P.append(dp_P)
    #order from method
    order_P.append(Order_P)
#Plot the stabilization Diagram    
plot = StabDia(nat_freqs_P, fRF_G, _, Freq, order_P, 'PolyMAX', 'no')

#%%

def FreqSeg(FRF, Freq, seg):
    import matplotlib.pyplot as plt 
    # Create figure and subplot 
    plt.figure(figsize=(17, 8))
    plt.xlim(0,max(Freq))
    plt.grid()  
    plt.xlabel("Frequency[Hz]")
    plt.ylabel("FRF")
    FRF = [abs(frf) for frf in FRF]
    plt.semilogy(Freq, FRF, label="FRF")
    # multiple segments
    plt.vlines(seg, ymin=0, ymax=max(FRF), colors='purple', ls='--', lw=2, label='EMA Segments')
    #plt.axvspan(0.0, seg, alpha=0.5, color='orange')
    plt.title('Modal Analysis Frequency Segment')
    plt.show()
    return
#%%
plot = FreqSeg(fRF_G, Freq, [5.0, 10.3, 13.5, 21.8, 25.52, 30.07, 39.0, 48.0, 55.0])

