%reset -f
#%%
from DataPreprocessing import DataPrep
from Rational_Polynomial_Fraction_Method import RFPM
from StabilizationDiagram import StabDia
from collections import defaultdict
import math
#Data Preprocessing
Iters = [[1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 6, 7, 8, 9], [10, 11, 12, 13, 14, 15, 16, 17, 18], [19, 20, 21, 22, 23, 24, 25, 26, 27]]
tests = ["BR_AR", "DS_TLE", "DS_RLE", "DS_CTE"]
output = {"BR_AR":[], "DS_TLE":[], "DS_RLE":[], "DS_CTE":[]}
Nat_Freq = []
AMP_levels = [[0.4, 0.8, 1.2, 1.6, 2], [0.4, 1.2, 2, 0.4, 1.2, 2, 0.4, 1.2, 2], [0.4, 1.2, 2, 0.4, 1.2, 2, 0.4, 1.2, 2], [0.4, 1.2, 2, 0.4, 1.2, 2, 0.4, 1.2, 2]]
DMG_levels = [[0, 0, 0, 0, 0], [256, 256, 256, 620, 620, 620, 898, 898, 898], [898, 898, 898, 4400, 4400, 4400, 8360, 8360, 8360], [256, 256, 256, 620, 620, 620, 898, 898, 898]]
for test_series, iters, AMP_level, DMG_level in zip(tests,Iters,AMP_levels,DMG_levels):
    reps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    iters_iter = iters
    AMP_level_iter =  AMP_level
    counter = [i for i in range(len(iters))]
    DMG_level_iter = DMG_level
    for iters, amp_levels, dmg_levels, i in zip(iters_iter, AMP_level_iter, DMG_level_iter, counter): 
        frf, freq = DataPrep([iters], reps, test_series)
        ranges = (
    ((5, 9), 1),
    ((12, 14), 1),
    ((16, 16.5), 1),
    ((15, 19), 1),
    ((18, 22), 1),
    ((22, 24), 2),
    ((26, 30), 1),
    ((30.5, 31), 1),
    ((35, 37), 1),
    ((40.5, 44), 2),
    ((48.5, 54), 2),)
        num_ord = 6
        import numpy as np
        #Experiment Description
        output[test_series].append({"Exp": test_series+"_"+str(iters), "AMP_level": amp_levels, "DMG_level": dmg_levels, "Wn":[], 'S.D': []})
        
        #RFPM parameters
        for (min_freq, max_freq), n_mode in ranges:
            wns = []
            wns1 = []
            for sensor_name in list(frf.keys()):
                #RFPM
                wn,_,_  = RFPM(frf, freq, min_freq, max_freq, sensor_name, n_mode, num_ord)
                #Natural frequency
                if n_mode ==2:
                    rounded_numbers = [num//1 for num in wn]
                    num_dict = defaultdict(list)
                    for num in wn:
                        num_dict[num//1].append(num)
                    modes = sorted(num_dict.values(), key=len, reverse=True)[:n_mode]
                    modes = sorted(modes, key=lambda x: x[0])
                    if len(modes) == 2:
                        wns.append(np.mean(modes[0]))
                        wns1.append(np.mean(modes[1]))
                    else:
                        wns.append(np.mean(modes[0][0]))
                        wns1.append(np.mean(modes[0][1]))
                else:
                    wns.append(np.mean(wn))
            if n_mode ==2:
                wns = list(filter(lambda x: not math.isnan(x), wns))
                wns1 = list(filter(lambda x: not math.isnan(x), wns1))
                output[test_series][i]["Wn"].append(round(np.mean(wns),3))
                output[test_series][i]["Wn"].append(round(np.mean(wns1),3))
                output[test_series][i]["S.D"].append(round(np.std(wns)*100,1))
                output[test_series][i]["S.D"].append(round(np.std(wns1)*100,1))
                #print(wns)
                #print(wns1)
            else:
                wns = list(filter(lambda x: not math.isnan(x), wns))
                output[test_series][i]["Wn"].append(round(np.mean(wns),3))
                output[test_series][i]["S.D"].append(round(np.std(wns)*100,1))
                #print(wns)
        Nat_Freq.append(output[test_series][i]["Wn"])
        Nat_Freq.append(output[test_series][i]["S.D"])
#%%

for ts in output.keys():
    
    for i in range(len(output[ts])):
        if ts == "DS_RLE" or ts== "DS_CTE":
            iters = [int(output[ts][i]["Exp"][-2:])]
        else:
            iters = [int(output[ts][i]["Exp"][-1:])]
        reps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        frf, freq = DataPrep(iters, reps, ts)
        frf_array = np.array(list(frf.values()))
        frf_mean = np.mean(frf_array, axis = 0)
        plot = FreqSeg(frf_mean, freq, output[ts][i]["Wn"], output[ts][i]["Exp"])

#%%
# Convert dictionary to LaTeX table format
for i in output.values():
    for j in i:
        latex_table = '\\begin{tabular}{|c|c|c|}\n\\hline\n'
        latex_table += 'Wn & S.D\\\\\n\\hline\n'
        for i in range(len(j['Wn'])):
            latex_table += f"{j['Wn'][i]} & {j['S.D'][i]}\\\\\n"
        latex_table += '\\hline\n\\end{tabular}'
        print(latex_table)
#%%
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Assume X is your data matrix

pca = PCA()
pca.fit(np.array(Nat_Freq))

# Plot the scree plot
plt.figure(figsize=(8, 8))
plt.plot(range(1, pca.n_components_+1), pca.explained_variance_ratio_, 'ro-')
plt.title('Scree Plot')
plt.xlabel('Principal Component')
plt.ylabel('Explained Variance Ratio')
plt.show()

# Determine the optimal number of components
cumulative_variance_ratio = np.cumsum(pca.explained_variance_ratio_)
num_components = np.argmax(cumulative_variance_ratio >= 0.95) + 1

pca = PCA(n_components=num_components)
pca.fit(np.array(Nat_Freq))
wn_pca = pca.transform(np.array(Nat_Freq))

# Create a list of markers that is repeated for each unique value of cat2
cat1_values = np.unique(sum(AMP_levels, []))
cat1_markers = ['o', 's', 'D', '^', 'v']  # Customize the marker shapes as needed
cat1_markers = cat1_markers[:len(sum(AMP_levels, []))] * (len(sum(AMP_levels, [])) // len(cat1_markers) + 1)


# Create a list of markers that is repeated for each unique value of cat2
cat2_values = np.unique(sum(DMG_levels, []))
cat2_colors = ['red', 'blue', 'green', 'orange', 'purple']  # Customize the marker colors as needed
cat2_colors = cat2_colors[:len(sum(DMG_levels, []))] * (len(sum(DMG_levels, [])) // len(cat2_colors) + 1)

# Plot the data points with different marker shapes and colors based on cat1 and cat2
for i in range(32):
    marker = cat1_markers[i]
    color = cat2_colors[i]
    plt.scatter(wn_pca[i, 0], wn_pca[i, 1], marker=marker, color=color)

# Add legend and axis labels
#plt.legend(handles=[plt.scatter([],[],marker=sum(AMP_levels, [])[i], color=sum(DMG_levels, [])[i], label=sum(AMP_levels, [])[i]) for i in range(32)])
#plt.legend(handles=[plt.scatter([],[],marker=cat1_dict[cat1_values[i]], color=cat1_color_dict[cat1_values[i]], label=cat1_values[i]) for i in range(len(cat1_values))])
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('First Two Principal Components')

# Show the plot
plt.show()
#%%

def FreqSeg(FRF, Freq, seg, test_series):
    import matplotlib.pyplot as plt 
    import numpy as np
    import pandas as pd
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
    plt.savefig('Results\\RFPM\\Out-of-band Estimate\\'+test_series+".png")
    plt.show()
    return
#%%
counter = 0 
for i in range(0, 59, 2):
    plot = FreqSeg(frf, freq, [5, 9, 12, 14, 22, 24, 26, 30, 30.5, 31, 35,37,40.5, 44, 48.5, 54],test_series,i, i+3, counter)
    counter +=1

