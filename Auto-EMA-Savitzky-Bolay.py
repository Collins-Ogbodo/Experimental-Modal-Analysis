from DataPreprocessing import DataPrep
from Rational_Polynomial_Fraction_Method import RFPM
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.signal import savgol_filter
#Data Preprocessing
iters = [27]
reps = [1,2,3,4,5,6,7,8,9,10]
test_series = "DS_CTE"
frf, freq, coh = DataPrep(iters, reps, test_series)
# Generate noisy data
min_freq, max_freq = 5, 60
FRF_df = pd.DataFrame(frf)
FRF_df_abs = abs(FRF_df)
imin = freq.index(min_freq)
imax = freq.index(max_freq)
FRF_df_abs = FRF_df_abs.iloc[imin:imax, :]
x = freq[imin:imax]
y = np.mean(FRF_df_abs, axis=1)
# Apply Savitzky-Golay smoothing
y_smooth = savgol_filter(y, window_length=35, polyorder=10)

plt.figure(figsize=(17, 8))
plt.xlim(min_freq,max(x))
plt.grid()  
plt.xlabel("Frequency[Hz]")
plt.ylabel("FRF")
plt.plot(x, y, 'b', label='Original')
plt.plot(x, y_smooth, 'r', label='Smoothed')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=10)
#plt.savefig('Results\\RFPM\\Out-of-band Estimate\\'+test_series+".png")
plt.show()

from scipy.signal import find_peaks, peak_prominences
peaks, _ = find_peaks(y_smooth, prominence=(1e-04))
prominences = peak_prominences(y_smooth, peaks)[0]
contour_heights = y_smooth[peaks] - prominences
plt.figure(figsize=(17, 8))
plt.grid() 
plt.plot( y_smooth)
plt.plot(peaks,  y_smooth[peaks], "x")
plt.vlines( x=peaks, ymin=contour_heights, ymax= y_smooth[peaks])
plt.show()
Id_peaks = [x[i] for i in peaks]
output = {}
incr =freq[1]
WNs = []

for i in Id_peaks:
    label = 'Mode'+str(Id_peaks.index(i)+1)
    label1 = 'STD-Mode'+str(Id_peaks.index(i)+1)
    output[label] = []
    output[label1] = []
    Wn = []
    for sensor_name in list(frf.keys()):
        freq_ini = [0, i]
        mmin_freq = i - incr
        mmax_freq = i + incr
        errors = {}
        
        for j in range(1,150):
            wn,_,_,_,_,_,error  = RFPM(frf, freq, mmin_freq, mmax_freq, sensor_name, n_modes=6)
            wn = np.mean(wn)
            if np.isnan(wn):
                wn = 0
            errors[str(wn)]=error
            mmin_freq-= incr
            mmax_freq+= incr
        best_wn = list(errors.keys())[list(errors.values()).index(min(errors.values()))]
        Wn.append(float(best_wn))
    output[label].append(np.mean(Wn))
    WNs.append(np.mean(Wn))
    output[label1].append(np.std(Wn))

def FreqSeg(FRF, Freq, seg, test_series):
    import matplotlib.pyplot as plt 
    import numpy as np
    import pandas as pd
    FRF = np.mean(np.array(list(FRF.values())), axis=(0))
    # Create figure and subplot
    imin = Freq.index(0.0)
    imax = Freq.index(60)
    #converting all data to array
    Freq = np.array(Freq)
    #Frequency, Coh and FRF range
    Freq = Freq[imin:imax]
    FRF = FRF[imin:imax]
    import numpy as np
    plt.figure(figsize=(17, 8))
    plt.xlim(0,max(Freq))
    plt.grid()  
    plt.xlabel("Frequency[Hz]")
    plt.ylabel("FRF")
    FRFs = [abs(frf) for frf in FRF]
    plt.semilogy(Freq, FRFs,color = 'orange', label = "Actual FRF")
    # multiple segments
    plt.vlines(seg, ymin=0, ymax=1,colors='purple', ls='--', lw=2, label='Estimate Natural Frequency')
    #for i in range(len(seg)):
       # plt.text(seg[i],seg[i+1], 'Segment'+str(i),  ha='center', va='bottom' )
    #plt.title('Modal Analysis RF-'+ test_series)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=10)
    #plt.savefig('Results\\RFPM\\Out-of-band Estimate\\'+test_series+".png")
    plt.show()
    return
plot_est = FreqSeg(frf, freq,WNs,test_series)

