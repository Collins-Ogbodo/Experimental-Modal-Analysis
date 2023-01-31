def PolyMAX(FRF, Freq, min_freq, max_freq, N):
    import numpy as np
    import itertools
    from scipy import signal
    import pandas as pd
    """This function computes the modal parameter of a system 
    using using the PolyMAX method FRF- Response Frequency as a  
    dictionary with keys as sensor name and values as list of FRF
    Frequency value as list Sensor input must be a the string of 
    the sensor name"""
    #N - degree of freedom
    #convert dictionary FRF to DataFrame
    FRF = pd.DataFrame(FRF)
    #convert dataframe to numpy arrary for splitting ease.
    FRF.to_numpy()
    #Extract the number of output
    m = FRF.shape[1]
    #Creating the FRF Matrix
    FRF = [list(frfs) for frfs in FRF.itertuples(index=False, name=None)]
    H = FRF
    print(FRF)
    #Because the frequency for all sensors is same whe just select the first.
    Freq = Freq['EXH']
    # Find corresponding indices of frequency range
    imin = Freq.index(min_freq)
    imax = Freq.index(max_freq) 
    #Frequency and FRF range
    Freq = Freq[imin:imax]
    FRF = FRF[imin:imax]
    #Length of frequency
    Nf = len(Freq)
    #compute the sample time
    f_0 = min(Freq); f_N = max(Freq)
    dt = 1/(2*(f_N-f_0))
    X_sensor = []
    for i in range(Nf):
        Omega = []
        for j in range(0,N+1):
            omega = np.exp(1j * Freq[i] * dt * j)
            Omega.append(omega)
        W_k = 1/Freq[i] if Freq[i] != 0.0 else 0.0
        X_ = np.multiply(W_k,Omega)
        X_sensor.append(X_)
    Y_sensor = np.kron(-X_,FRF)
    R_sensor = np.real(np.dot(np.transpose(np.conjugate(X_sensor)), X_sensor))
    T_sensor = np.real(np.dot(np.transpose(np.conjugate(Y_sensor)), Y_sensor))
    S_sensor = np.real(np.dot(np.transpose(np.conjugate(X_sensor)), Y_sensor)) 
    M_sensor = T_sensor - (np.transpose(np.conjugate(S_sensor)) @ np.linalg.inv(R_sensor) @ S_sensor)
    M_sensor = 2* M_sensor
    alpha = -1*np.linalg.inv(M[0:N*m,0:N*m]) @ M[0:N*m, N*m:len(M)]
    Im = np.identity(m)
    alpha = np.concatenate((alpha, Im), axis = 0)
    #creatin the companion matrix
    CM = np.eye(len(M),k = 1)
    #remove the last row and join with alpha
    CM = CM[0:-1,:]
    CM = np.concatenate((CM, np.transpose(alpha)), axis = 0)
    #poles, PF = np.linalg.eig(CM)
    #poles = np.log(poles)/dt
    #stable poles
    #poles = [pole for pole in poles if np.real(pole) < 0.0 and np.imag(pole) > 0.0]
    #Compuring Natural Frequency
    #w_n = [abs(pole) for pole in poles]
    return M_sensor, alpha, CM,H
    
#%%
from DataPreprocessing import DataPrep
#Data Preprocessing
iters = [1]
reps = [1]
test_series = "BR_AR"
frf, freq = DataPrep(iters, reps, test_series)

M, alpha, cm, H = PolyMAX(frf, freq, 47, 55, 20)





