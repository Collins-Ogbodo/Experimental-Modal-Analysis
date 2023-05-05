from DataPreprocessing import DataPrep
from Rational_Polynomial_Fraction_Method import RFPM
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
#Data Preprocessing
iters = [27]
reps = [1]
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
peaks_all, _ = find_peaks(y_smooth, prominence=(0))
prominences_all = peak_prominences(y_smooth, peaks_all)[0]
contour_heights_all = y_smooth[peaks_all] - prominences_all
plt.figure(figsize=(17, 8))
plt.grid() 
plt.title("All Identified Peaks")
plt.plot( y_smooth)
plt.plot(peaks_all,  y_smooth[peaks_all], "x")
plt.vlines( x=peaks_all, ymin=contour_heights_all, ymax= y_smooth[peaks_all])
plt.show()

peaks_auto, _ = find_peaks(y_smooth, prominence=(1e-04))
prominences_auto = peak_prominences(y_smooth, peaks_auto)[0]
contour_heights_auto = y_smooth[peaks_auto] - prominences_auto
plt.figure(figsize=(20, 8))
plt.grid() 
plt.title("Auto-Selected Peaks")
plt.plot( y_smooth)
plt.plot(peaks_auto,  y_smooth[peaks_auto], "x")
plt.vlines( x=peaks_auto, ymin=contour_heights_auto, ymax= y_smooth[peaks_auto])
plt.show()

peaks_add, _ = find_peaks(y_smooth, prominence=(0))
prominences_add = peak_prominences(y_smooth, peaks_add)[0]
contour_heights_add = y_smooth[peaks_add] - prominences_add
# Define onclick function
def onclick(event, peak):
    # Get the x-coordinate of the click
    p = event.xdata
    if p is not None:
        # Find the index of the closest point in the data array
        idx = (np.abs(p - peaks_add)).argmin()
        print(idx)
        # Add the x-coordinate and y-value to the peak list
        peak.append(peaks_add[idx])
        print("Selected peak:", peaks_add[idx])
 
     
fig, ax = plt.subplots(figsize=(17, 8))
ax.grid() 
ax.plot( y_smooth)
plt.title("Select New Peaks Here")
ax.plot(peaks_add,  y_smooth[peaks_add], "x")
ax.vlines( x=peaks_add, ymin=contour_heights_add, ymax= y_smooth[peaks_add])
# Create list to hold selected peaks
Id_peaks_add = []
# Connect onclick function to plot
cid = fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, Id_peaks_add))
plt.show()


def FreqSeg(FRF, Freq, seg, test_series, min_freq, max_freq):
    import matplotlib.pyplot as plt 
    import numpy as np
    FRF = np.mean(np.array(list(FRF.values())), axis=(0))
    # Create figure and subplot
    imin = Freq.index(min_freq)
    imax = Freq.index(max_freq)
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
#%%
f_peaks = np.unique(list(peaks_auto) + list(Id_peaks_add))
Id_peaks = [ x[i] for i in f_peaks]
prominences = peak_prominences(y_smooth, f_peaks)[0]
contour_heights = y_smooth[f_peaks] - prominences
plt.figure(figsize=(17, 8))
plt.title("Final Peaks Idenfied")
plt.grid()
plt.plot( y_smooth)
plt.plot(f_peaks,  y_smooth[f_peaks], "x")
plt.vlines( x=f_peaks, ymin=contour_heights, ymax= y_smooth[f_peaks])
plt.show()
#%%
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
        errors_o ={}
        for k in range(5,20):
            freq_ini = [0, i]
            mmin_freq = i - incr
            mmax_freq = i + incr
            errors = {}
            for j in range(1,150):
                wn,_,_,_,_,_,error  = RFPM(frf, freq, mmin_freq, mmax_freq, sensor_name, n_modes=k)
                wn = np.mean(wn)
                if np.isnan(wn):
                    wn = 0
                errors[str(wn)]=error
                mmin_freq-= incr
                mmax_freq+= incr
            min_error = min(errors.values())
            best_wn = list(errors.keys())[list(errors.values()).index(min_error)]
            errors_o[min_error] = float(best_wn)
        best_wn_o = errors_o[min(errors_o.keys())]
        Wn.append(float(best_wn_o))
    output[label].append(np.mean(Wn))
    WNs.append(np.mean(Wn))
    output[label1].append(np.std(Wn))

plot_est = FreqSeg(frf, freq,WNs,test_series,  min_freq, max_freq)

