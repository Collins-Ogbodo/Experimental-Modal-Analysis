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
ranges =((5, 9),
    (12, 14),
    (15, 16.5),
    (15, 19),
    (18, 22),
    (22, 24),
    (26, 30),
    (30.5, 31),
    (35, 37),
    (40.5, 44),
    (48.5, 54),)
for i, j in ranges:
    # Generate noisy data
    min_freq, max_freq = i, j
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
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=10)
    #plt.savefig('Results\\RFPM\\Out-of-band Estimate\\'+test_series+".png")
    plt.show()
