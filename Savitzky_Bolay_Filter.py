import pandas as pd
from scipy.signal import savgol_filter
import numpy as np
import matplotlib.pyplot as plt
FRF_filter = pd.DataFrame(frf)
FRF_filter = FRF_filter.to_numpy()

actual_FRF = FRF_filter
order = [2, 5, 9]
window = [3, 6, 9]
for o in order:
    window = np.array(window)+o
    for w in window:    
        FRF_filter = savgol_filter(FRF_filter, w, o, axis = 1, mode="interp" )
        plt.figure(figsize=(17, 8))
        plt.xlim(0,max(freq[0:1920]))
        plt.semilogy(freq[0:1920],abs(FRF_filter[0:1920,1]), label = 'Filtered')
        plt.semilogy(freq[0:1920],abs(actual_FRF[0:1920,1]), label = 'actual')
        
        plt.legend()
        plt.text(58, 0.006, "Order-"+str(o), horizontalalignment='center',
         verticalalignment='center')
        plt.text(58, 0.009, "Window-"+str(w), horizontalalignment='center',
         verticalalignment='center')
        plt.title(test_series+"-EXH")
        plt.show()