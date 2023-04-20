from DataPreprocessing import DataPrep
from Rational_Polynomial_Fraction_Method import RFPM
import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
#Data Preprocessing
iters = [1]
reps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
test_series = "BR_AR"
min_freq, max_freq = 0, 60
frf, freq, coh = DataPrep(iters, reps, test_series)
#%%
FRF_df = pd.DataFrame(frf)
FRF_df = abs(FRF_df)
imin = freq.index(min_freq)
imax = freq.index(max_freq)
FRF_df = FRF_df.iloc[imin:imax, :]
#%%
FRF_df.info()
FRF_df.isnull
#%%
kmeans = KMeans(n_clusters= 14, init='k-means++', random_state=0)
kmeans.fit(FRF_df)
label = kmeans.labels_
wcss = kmeans.inertia_
    
#%%

plt.scatter(FRF_df.iloc[:,50],FRF_df.iloc[:,2], c= kmeans.labels_)
