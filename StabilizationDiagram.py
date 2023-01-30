def StabDia(NatFreq, FRF, FRF_est, Freq, Order, sensor):
    import matplotlib.pyplot as plt 
    # Create figure and subplot
    fig, host = plt.subplots(figsize=(15,10))   
    plt.grid()  
    par1 = host.twinx()
    host.set_xlabel("Frequency[Hz]")
    host.set_ylabel("FRF")
    par1.set_ylabel("Model Order")
    par1.set_yticks(Order)
    color1 = plt.cm.viridis(0)
    FRF = [abs(frf) for frf in FRF]
    host.semilogy(Freq, FRF, color=color1, label="FRF")
    for i, j in zip(FRF_est, Order):
        k = [abs(frf) for frf in i]
        host.semilogy(Freq, k)
    host.legend(Order, loc="upper left")
    for w_n, N in zip(NatFreq,Order):
        N= [N for i in range(len(w_n))]
        par1.plot(w_n, N,"+", markersize=5, label='gaga')
    plt.title(sensor)
    plt.show()
    return fig