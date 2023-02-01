def GRFPM(FRF, Freq, min_freq, max_freq, N):
    import numpy as np
    import itertools
    from scipy import signal
    import pandas as pd
    """This function computes the modal parameter of a system 
    using the Rational Fraction Polynomial Method
    FRF- Response Frequency as a dictionary with keys as sensor 
    name and values as list of FRF
    Frequency value as list 
    Sensor input must be a the string of the sensor name"""
    #N - degree of freedom
    n = 2*N  #number of denominator polynomial terms
    #convert dictionary FRF to DataFrame
    FRF = pd.DataFrame(FRF)
    #find the average of all sensors append this average in the last column
    FRF= FRF.mean(axis=1)
    #convert dataframe to numpy arrary for splitting ease.
    FRF.to_numpy()
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
    V = np.linalg.multi_dot([np.transpose(X),np.linalg.inv(Y),X]) + Z
    U = np.linalg.multi_dot([np.transpose(X),np.linalg.inv(Y),G]) - F
    #obtaining the coefficient of of the characteristics polynomial
    b = np.dot(np.linalg.inv(V),U)
    #calculating the coefficient of the numerator of the tranfer function.
    a = np.dot(np.linalg.inv(Y),(G-np.dot(X,b)))
    #convert the dimension of a and b to vectors
    a = list(itertools.chain.from_iterable(a))
    b = list(itertools.chain.from_iterable(b))
    #Reverse the order of coefficient
    a.reverse(), b.reverse()
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


