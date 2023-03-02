def StabDia(NatFreq, FRF, FRF_est, Freq, Order, sensor, recon = 'yes', algo =''):
    import matplotlib.pyplot as plt 
    import numpy as np
    # Create figure and subplot
    fig, host = plt.subplots(figsize=(15,10))   
    plt.grid()  
    par1 = host.twinx()
    host.set_xlabel("Frequency[Hz]")
    host.set_ylabel("FRF")
    par1.set_ylabel("Model Order")
    par1.set_yticks(Order)
    color1 = plt.cm.viridis(0)
    if algo =='P':
        for i in range(np.shape(FRF)[2]):
            FRFs = [abs(frf) for frf in FRF[0,:,i]]
            host.semilogy(Freq, FRFs)
    else:
            FRF = [abs(frf) for frf in FRF]
            host.semilogy(Freq, FRF, color=color1, label="FRF")
        
    if recon == 'yes':
        for i, j in zip(FRF_est, Order):
            k = [abs(frf) for frf in i]
            host.semilogy(Freq, k, label =j)
    elif recon == 'no':
        pass
    else:
        print('recon = yes/no')
    host.legend(loc="upper left")
    for w_n, N in zip(NatFreq,Order):
        N= [N for i in range(len(w_n))]
        #par1.plot(w_n, N,"+", markersize=5)
        if len(w_n) > 0:
            plt.vlines(w_n, ymin=0, ymax=1, ls='--', lw=2, label='wn')        
    plt.title(sensor)
    plt.show()
    return fig