
def RFPM(FRF, Freq, min_freq, max_freq, sensor, n_modes, num_ord):
    import numpy as np
    """This function computes the modal parameter of a system 
    using the Rational Fraction Polynomial Method
    FRF- Response Frequency as a dictionary with keys as sensor 
    name and values as list of FRF
    Frequency value as list 
    Sensor input must be a the string of the sensor name"""
    #Number of coefficient in the Numerator
    m = n_modes * 2 + 1 + num_ord
    #number of denominator polynomial terms
    n = n_modes * 2 + 1
    #Selecting Sensor
    FRF = FRF[sensor]
    # Find corresponding indices of frequency range
    imin = Freq.index(min_freq)
    imax = Freq.index(max_freq) 
    #Frequency and FRF range
    w = np.array(Freq[imin:imax])
    H = np.array(FRF[imin:imax])
    # Hermitian transpose
    def HT(a):
        return a.conj().T
    # complex frequency vector
    iw = 1j * w

    # Build monomial basis matricies
    Phi_a = iw[:, None] ** np.arange(m)
    Phi_b_all = iw[:, None] ** np.arange(n)
    Phi_b = Phi_b_all[:, :-1]  # ignore last column because bn=1

    # preallocate some calculations for speed
    Phi_bH = Phi_b * H[:, None]
    Hiwn = H * Phi_b_all[:, -1]
    D = -HT(Phi_a) @ (Phi_bH)

    # form the block matricies
    M = np.block([[HT(Phi_a) @ Phi_a, D], [D.T, HT(Phi_bH) @ (Phi_bH)]])
    x = np.block([HT(Phi_a) @ Hiwn, -HT(Phi_bH) @ Hiwn])

    # Solve and extract the coefficients of the polynomials
    AB = np.linalg.solve(np.real(M), np.real(x))
    a = AB[:m, None]
    b = np.append(AB[m:], 1)[:, None]

    # Generate the predicted FRF
    H_pred = (Phi_a @ a) / (Phi_b_all @ b)

    # Pull out the modal porperties
    roots_b = sorted(np.roots(np.flip(b[:, 0])))[
        ::-2
    ]  # remove every other becaus they are conj pairs
    wns = np.abs(roots_b)
    zetas = -np.real(roots_b) / wns
    return wns, H_pred, zetas

