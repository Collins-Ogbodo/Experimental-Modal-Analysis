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
    #Because the frequency for all sensors is same whe just select the first.
    Freq = Freq['EXH']
    # Find corresponding indices of frequency range
    imin = Freq.index(min_freq)
    imax = Freq.index(max_freq) 
    #Frequency and FRF range
    Freq = Freq[imin:imax]
    #Length of frequency
    Nf = len(Freq)
    #compute the sample time
    f_0 = min(Freq); f_N = max(Freq)
    #Sampling Time
    dt = 1/(2*(f_N-f_0))
    #Define initial matrix for M
    M = np.zeros([N+1,N+1])
    #Number ouput
    l = len(FRF)
    #number input
    m = len(FRF[0])
    X = []
    for i in range(Nf):
        Omega = []
        for j in range(0,N+1):
            omega = np.exp(1j * Freq[i] * dt * j)
            Omega.append(omega)
        #print(Freq[i])
        W_k = 1/Freq[i] if Freq[i] != 0.0 else 0.0
        X_sensor = np.multiply(W_k,Omega)
        X.append(X_sensor)
    #compute negative X
    #X_neg = [[-1*x for x in innerlist] for innerlist in X]
    #compute M for all output
    for ls in range(l):
        #Repeat for number of input and compute the Y columns
        Y = [] 
        for ms in range(m):
            H_sensor = FRF[ls][ms]
            H_sensor = H_sensor[imin:imax]
            for x, h in zip(X,H_sensor):
                X_neg = [-1*xs for xs in x]
                Y_x = np.multiply(X_neg,h)
                Y_x = list(Y_x)
                Y.append(Y_x)
        R_sensor = np.real(np.dot(np.transpose(np.conjugate(X)), X))
        T_sensor = np.real(np.dot(np.transpose(np.conjugate(Y)), Y))
        S_sensor = np.real(np.dot(np.transpose(np.conjugate(X)), Y)) 
        M_sensor = T_sensor - (np.transpose(np.conjugate(S_sensor)) @ np.linalg.inv(R_sensor) @ S_sensor)
        M =+ M_sensor
    M = 2*M
    alpha = -1*np.linalg.inv(M[0:N*m,0:N*m]) @ M[0:N*m, N*m:len(M)]
    Im = np.identity(m)
    alpha = np.concatenate((alpha, Im), axis = 0)
    #creating the companion matrix
    CM = np.eye(len(M),k = 1)
    #remove the last row and join with alpha
    CM = CM[0:-1,:]
    CM = np.concatenate((CM, np.transpose(alpha)), axis = 0)
    poles, PF = np.linalg.eig(CM)
    poles = np.log(poles)/dt
    stable poles
    #poles = [pole for pole in poles if np.real(pole) < 0.0 and np.imag(pole) > 0.0]
    #Compuring Natural Frequency
    #w_n = [abs(pole) for pole in poles]
    return M, alpha, CM, poles
  
#%%


def PolyMaxDataPrep(frf):
    """This function is pecific to the data sent available, it coverts the two dimensional 
    space data (SIMO) which is a dictionary with keys equal to the sensor (output)
    and vales equal to a list of each input across all frequencies to a three dimensional
    matix space
    e.g [[[FRF(w) for input 1 for ouput 1 ], [FRF(w) for input 2 for ouput 1 ], ...],
         [[FRF(w) for input 1 for ouput 2 ], [FRF(w) for input 2 for ouput 2 ], ...]], ...]"""
    FRF = []
    for i in frf.values():
        dump_dim = []
        dump_dim.append(i)
        FRF.append(dump_dim)
    return FRF
          
#%%
from DataPreprocessing import DataPrep
#Data Preprocessing
iters = [1]
reps = [1]
test_series = "BR_AR"
frf, freq = DataPrep(iters, reps, test_series)

fRF = PolyMaxDataPrep(frf)


#%%

M, alpha, cM, p = PolyMAX(fRF, freq, 47, 55, 20)





