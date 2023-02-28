import pandas as pd
from scipy.signal import savgol_filter
import numpy as np
import matplotlib.pyplot as plt
FRF_filter = pd.DataFrame(frf)
FRF_filter = FRF_filter.to_numpy()

actual_FRF = FRF_filter
frf_real = np.real(FRF_filter)
frf_imag = np.imag(FRF_filter)
FRF_filter_real = savgol_filter(frf_real, 10, 2, axis = 1, mode="nearest" )
FRF_filter_imag = savgol_filter(frf_imag, 10, 2, axis = 1, mode="nearest" )
FRF_filter = FRF_filter_real +1j * FRF_filter_imag

plt.semilogy(freq[0:500],abs(FRF_filter[0:500,1]), label = 'Filtered')
plt.semilogy(freq[0:500],abs(actual_FRF[0:500,1]), label = 'actual')
plt.legend()
plt.show()