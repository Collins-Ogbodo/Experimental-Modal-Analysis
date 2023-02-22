def PolyMAX(FRF, Freq, Coh, min_freq, max_freq, Nmin, Nmax):
    import numpy as np
    """This function computes the modal parameter of a system 
    using using the PolyMAX method. The FRF input is a 3 dimensional matrix
    of the output, input and corresponding frequency
    Frequency value is a dictionary of a frequencies or output"""
    #N - degree of freedom
    N = list(range(Nmax, Nmin-1,-1))
    #Find the maximum Order
    Nmax = max(N)
    # Find corresponding indices of frequency range
    imin = Freq.index(min_freq)
    imax = Freq.index(max_freq)
    #converting all data to array
    Freq = np.array(Freq)
    #Frequency, Coh and FRF range
    Freq = Freq[imin:imax]
    FRF = FRF[:, imin:imax, :]
    #Coh = Coh[:, imin:imax, :]
    #Length of frequency-Nf, number of input-m, number of output-l
    m, Nf, l = np.shape(FRF)
    #compute the sample time
    f_0 = min(Freq); f_N = max(Freq)
    #Sampling Time
    dt = 1/(2*(f_N-f_0))
    #Define initial matrix for M, X, Y, for maximum model order 0 to N
    M = np.zeros([Nmax+1,Nmax+1], dtype=complex)
    X = np.zeros([Nf, Nmax+1], dtype=complex)
    Y = np.zeros([Nf, (Nmax+1)*m], dtype=complex)
    #Scalar weighting function
    W_k = 1#np.array([1/freq_val if freq_val != 0.0 else 0.0 for freq_val in Freq])
    #W_k = 100 #Coh[ms][i][ls]
    #computing X
    for j in range(Nf):
        x = np.exp(1j * Freq[j] * dt * np.array([range(Nmax+1)])) * W_k#[j]
        X[j,:] = x 
    R_sensor = np.real(np.matmul(np.transpose(X), X))
    #compute M for all output
    for ls in range(l):
        for i in range(Nf):
            Y[i,:] = -X[i,:] * np.transpose(FRF[:,i,ls])
        T_sensor = np.real(np.matmul(np.transpose(Y), Y))
        S_sensor = np.real(np.matmul(np.transpose(X), Y)) 
        M_sensor = T_sensor - (np.transpose(S_sensor) @ np.linalg.inv(R_sensor) @ S_sensor)
        M = M + M_sensor 
    M = 2*M
    #Iterating for all model order
    nat_freqs_P =[]
    damp_ratio_P = []
    order_P = []
    for n in N:
        A = M[0:n*m,0:n*m]
        B = M[0:n*m, (n*m):]
        alpha,_,_,_ = np.linalg.lstsq(-A, np.reshape(B, (-1,1)), rcond=None)
        Im = np.identity(m)
        alpha = np.concatenate((alpha, Im), axis = 0)
        #creating the companion matrix
        CM = np.eye(len(M),k = 1)
        #remove the last row and join with alpha
        CM = CM[0:-1,:]
        CM = np.concatenate((CM, -np.transpose(alpha)), axis = 0)
        poles, PF = np.linalg.eig(CM)
        poles = np.log(poles)/dt
        # Picking only poles with stable(negative) poles
        poles = [pole for pole in poles if np.real(pole) < 0.0 and np.imag(pole) > 0.0 
                 and np.abs(pole) <= max_freq and np.abs(pole) >= min_freq ]
        # Find the natural frequency and damping factor of the system
        nat_freq = np.abs(poles)
        #Computing the Damping Ratio
        dam_ratio = -np.real(poles) / nat_freq 
        #Natural frequency
        nat_freqs_P.append(nat_freq)  
        #Damping Ratio
        damp_ratio_P.append(dam_ratio)
        #order from method
        order_P.append(n)  
        M = M[:-1,:-1]
    return nat_freqs_P, damp_ratio_P,order_P





