def PolyMAX(FRF, Freq, min_freq, max_freq, N):
    import numpy as np
    """This function computes the modal parameter of a system 
    using using the PolyMAX method. The FRF input is a 3 dimensional matrix
    of the output, input and corresponding frequency
    Frequency value is a dictionary of a frequencies or output"""
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
    #Number output
    l = len(FRF)
    #number input
    m = len(FRF[0])
    X = []
    for i in range(Nf):
        Omega = []
        for j in range(0,N+1):
            omega = np.exp(-1j * Freq[i] * dt * j)
            Omega.append(omega)
        #print(Freq[i])
        W_k = 1/Freq[i] if Freq[i] != 0.0 else 0.0
        X_sensor = np.multiply(W_k,Omega)
        X.append(X_sensor)
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
    poles = -1*np.log(poles)/dt
    print(poles)
    # Picking only poles with stable(negative) poles
    poles = [pole for pole in poles if np.real(pole) < 0.0 and np.imag(pole) > 0.0 
             and np.abs(pole) <= max_freq and np.abs(pole) >= min_freq ]
    # Find the natural frequency and damping factor of the system
    nat_freq = np.abs(poles)
    #Computing the Damping Ratio
    dam_ratio = -np.real(poles) / nat_freq 
    return nat_freq, dam_ratio, N, FRF, Freq





