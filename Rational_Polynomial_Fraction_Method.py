def RFPM_alg(FRF, Freq, min_freq, max_freq, sensor, N):
    import numpy as np
    import itertools
    from scipy import signal
    """This function computes the modal parameter of a system 
    using the Rational Fraction Polynomial Method
    FRF- Response Frequency as a dictionary with keys as sensor 
    name and values as list of FRF
    Frequency value as list 
    Sensor input must be a the string of the sensor name"""
    #N - degree of freedom
    n = 2*N  #number of denominator polynomial terms
    #Selecting Sensor
    FRF = FRF[sensor]
    Freq = Freq[sensor]
    # Find corresponding indices of frequency range
    imin = Freq.index(min_freq)
    imax = Freq.index(max_freq) 
    #Frequency and FRF range
    Freq = Freq[imin:imax]
    FRF = FRF[imin:imax]
    #Generating individual Matrices
    P = []
    for i in Freq:
        new_row = []
        for k in range(n):
            new_ele= (complex(0,i))**k
            new_row.append(new_ele)
        P.append(new_row)
    T =[]
    for i, j in zip(Freq, FRF):
        new_row = []
        for k in range(n):
            new_ele= j*(complex(0,i))**k
            new_row.append(new_ele)
        T.append(new_row) 
    W= []
    for i, j in zip(Freq, FRF):
        new_row= [j*(complex(0,i))**n]
        W.append(new_row)
    Y = np.real(np.dot(np.transpose(np.conjugate(P)),P))
    X = -1*np.real(np.dot(np.transpose(np.conjugate(P)),T))
    Z = np.real(np.dot(np.transpose(np.conjugate(T)),T))
    G = np.real(np.dot(np.transpose(np.conjugate(P)),W))
    F = -1*np.real(np.dot(np.transpose(np.conjugate(T)),W))
    # Concatenation 
    top_half = np.concatenate((Y, X), axis = 1)
    but_half = np.concatenate((np.transpose(X), Z), axis = 1)
    matrix = np.concatenate((top_half, but_half), axis = 0)
    #comput inverse and find the right side vector
    vec = np.concatenate((G, F), axis = 0)
    coeff = np.dot(np.linalg.inv(matrix), vec)
    #Split coefficient into A and B
    a, b = np.array_split(coeff, 2)
    #convert the dimension of a and b to vectors
    a = list(itertools.chain.from_iterable(a)) 
    b = list(itertools.chain.from_iterable(b))
    #Reverse the order of the list to suit the residue function
    a.reverse(), b.reverse()
    #insert the coefficient of the variable  of highest power
    b = [1.0,*b]   
    # Solve the characteristics equation 
    poles = np.roots(b)
    # Picking only poles with stable(negative) poles
    poles = [pole for pole in poles if np.real(pole) < 0.0 and np.imag(pole) > 0.0 
             and np.abs(pole) <= max_freq and np.abs(pole) >= min_freq ]
    # Find the natural frequency and damping factor of the system
    nat_freq = np.abs(poles)
    #Computing the Damping Ratio
    dam_ratio = -np.real(poles) / nat_freq 
    #Create the FRF using the estimated coefficients
    _, FRF_est = signal.freqs(a, b, worN=Freq)
    return nat_freq, dam_ratio, n, FRF, Freq, FRF_est




