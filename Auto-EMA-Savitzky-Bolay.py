from DataPreprocessing import DataPrep
from Rational_Polynomial_Fraction_Method import RFPM
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.signal import savgol_filter
#Data Preprocessing
iters = [1]
reps = [1]
test_series = "BR_AR"
frf, freq, coh = DataPrep(iters, reps, test_series)
#%%
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
plt.plot( y_smooth)
plt.plot(peaks,  y_smooth[peaks], "x")
plt.vlines( x=peaks, ymin=contour_heights, ymax= y_smooth[peaks])
plt.show()
Id_peaks = [x[i] for i in peaks]

for i in Id_peaks:
    #Because the frquency spectrum is equally spaced we use this equal spacing
    incr = freq[1]
    mmin_freq = i - find_peaks
    mmax_freq = i + find_peaks
    
    #RFPM parameters
    wn,_,_,_,_,_  = RFPM(frf, freq, min_freq, max_freq, sensor_name, n_mode, num_ord)
