import numpy as np
import matplotlib.pyplot as plt
import itertools
import pandas as pd
m = 1
k =3.948e+5
w_n = [100, 50, 25]
damp_ratio = 0.004
Freq = [i for i in range(0,200)]
frf = []
for i in Freq:
    #frfs = (1/k)*((-1*(i**2)*(w_n**2))/((w_n**2)-(i**2)+1j*(2*damp_ratio*i*w_n)))
    frfs = (-1/m)*(pow(i,2)/np.sqrt((pow(w_n[0],2)-pow(i,2))**2+pow((2*damp_ratio*w_n[0]*i),2)))+\
        (-1/m)*(pow(i,2)/np.sqrt((pow(w_n[1],2)-pow(i,2))**2+pow((2*damp_ratio*w_n[1]*i),2)))+\
            (-1/m)*(pow(i,2)/np.sqrt((pow(w_n[2],2)-pow(i,2))**2+pow((2*damp_ratio*w_n[2]*i),2)))
    frf.append(frfs)
FREQ = {}
FRF = {}
FREQ['Freq'] = Freq
FRF['frf'] = frf

import scipy.io
scipy.io.savemat('Freq.mat', FREQ)
scipy.io.savemat('FRF.mat', FRF)
frf = pd.DataFrame(frf)
frf = frf.to_numpy()
def PolyMAX(frf, freq, N):
    """This function computes the modal parameters using the PolyMAX method"""
    
    #No. of sensors
    sens_no = frf.shape[1]
    
    #Length of frequency
    len_ = len(freq)

    #computing the base function for all frequency and model order
    #number of degree of freedom
    M = np.zeros([N+1,N+1])
    #compute the sample time
    f_0 = min(freq); f_N = max(freq)
    dt = 1/(2*(f_N-f_0))
    
    for sensor in range(sens_no):
        X_sensor = []
        Y_sensor = []
        for i in range(len_):
            Omega = []
            for j in range(0,N+1):
                omega = np.exp(1j * freq[i] * dt * j)
                Omega.append(omega)
            W_k = 1/freq[i] if freq[i] != 0.0 else 0.0
            X_ = np.multiply(W_k,Omega)
            X_sensor.append(X_)
            Y_sensor.append(np.multiply(-X_,frf[i,sensor]))
        R_sensor = np.real(np.dot(np.transpose(np.conjugate(X_sensor)), X_sensor))
        T_sensor = np.real(np.dot(np.transpose(np.conjugate(Y_sensor)), Y_sensor))
        S_sensor = np.real(np.dot(np.transpose(np.conjugate(X_sensor)), Y_sensor))
        M_sensor = T_sensor - (np.transpose(np.conjugate(S_sensor)) @ np.linalg.inv(R_sensor) @ S_sensor)
        M += M_sensor
    M = np.multiply(M,2)
    alpha = -1*np.linalg.inv(M[0:-1,0:-1]) @ [M[k,-1] for k in range(N)]
    alpha = np.append(alpha, 1.0)
    alpha = [k for k in alpha]
    #alpha = alpha[::-1]
    #creatin the companion matrix
    CM = np.eye(N,k = 1)
    CM[-1,:] = alpha[0:-1]
    poles, PF = np.linalg.eig(CM)
    poles = np.log(poles)/dt
    #stable poles
    poles = [pole for pole in poles if np.real(pole) < 0.0 and np.imag(pole) > 0.0]
    #Compuring Natural Frequency
    w_n = [abs(pole) for pole in poles]
    return w_n, M, CM, R_sensor, M, Omega
    
W_n, M, MC, R, M, Omega = PolyMAX(frf, Freq, 12)

#%%
