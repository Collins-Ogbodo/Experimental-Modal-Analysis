import numpy as np
import matplotlib.pyplot as plt
import itertools
m = 1
k =3.948e+5
w_n = [100, 50, 25]
damp_ratio = 0.004
Freq = [i for i in range(1,200)]
frf = []
for i in Freq:
    #frfs = (1/k)*((-1*(i**2)*(w_n**2))/((w_n**2)-(i**2)+1j*(2*damp_ratio*i*w_n)))
    frfs = (-1/m)*(pow(i,2)/np.sqrt((pow(w_n[0],2)-pow(i,2))**2+pow((2*damp_ratio*w_n[0]*i),2)))+\
        (-1/m)*(pow(i,2)/np.sqrt((pow(w_n[1],2)-pow(i,2))**2+pow((2*damp_ratio*w_n[1]*i),2)))+\
            (-1/m)*(pow(i,2)/np.sqrt((pow(w_n[2],2)-pow(i,2))**2+pow((2*damp_ratio*w_n[2]*i),2)))
    frf.append(frfs)
frf_ = [abs(i) for i in frf]

def G_RFPM_alg(FRF, Freq, N):
    """This function computes the modal parameter of a system 
    using the Rational Fraction Polynomial Method"""
    #N - degree of freedom
    
    n = 2*N    #number of denominator polynomial terms
    
    #Generating individual Matrices
    P = []
    for i in Freq:
        new_row = []
        for k in range(n):
            new_ele= (complex(0,i))**k
            new_row.append(new_ele)
        P.append(new_row)
    
    T =[]
    for i, j in zip(Freq, frf):
        new_row = []
        for k in range(n):
            new_ele= j*(complex(0,i))**k
            new_row.append(new_ele)
        T.append(new_row)
    
    W= []
    for i, j in zip(Freq, frf):
        new_row= [j*(complex(0,i))**n]
        W.append(new_row)
    
    Y = np.real(np.dot(np.transpose(np.conjugate(P)),P))
    X = -1*np.real(np.dot(np.transpose(np.conjugate(P)),T))
    Z = np.real(np.dot(np.transpose(np.conjugate(T)),T))
    G = np.real(np.dot(np.transpose(np.conjugate(P)),W))
    F = -1*np.real(np.dot(np.transpose(np.conjugate(T)),W))
    
    V = np.linalg.multi_dot([np.transpose(X),np.linalg.inv(Y),X]) + Z
    U = np.linalg.multi_dot([np.transpose(X),np.linalg.inv(Y),G]) - F
    
    coeff = np.dot(np.linalg.inv(V),U)
    
    b = coeff
    #convert the dimension of a and b to vectors
    b = list(itertools.chain.from_iterable(b))
    #Reverse the order of the list to suit the residue function
    b.reverse()
    b = [1.0,*b]
        
    # Solve the characteristics equation 
    poles = np.roots(b)
    
    # Picking only poles with stable(negative) poles
    poles = [p for p in poles if np.real(p) <= 0.0]
    
    # Find the natural frequency and damping factor of the system
    nat_freq = np.abs(poles)
    #Computing the Damping Ratio
    dam_ratio = -np.real(poles) / nat_freq
    
    #Remove frequencies more than the frequency spectrum under consideration
    nat_freq = [i for i in nat_freq]
    nat_freq = [x for x in nat_freq if x < max(Freq) and x >= min(Freq)]
    
    #remove damping factors outside the 0-1 range
    dam_ratio = [dr for dr in dam_ratio if dr < 1.0 and dr > 0.0]
    return nat_freq, dam_ratio, b, poles, n
N = [i for i in range(1,20)]
nat_freqs =[]
for i in N:
    nat_freq_ = RFPM_alg(frf, Freq, i)
    nat_freqs.append(nat_freq_)
    


frf_ = [abs(i) for i in frf]


import matplotlib.pyplot as plt 

# Create figure and subplot manually

fig, host = plt.subplots(figsize=(8,5))   
plt.grid()  
par1 = host.twinx()
    
host.set_xlabel("Frequency[Hz]")
host.set_ylabel("FRF")
par1.set_ylabel("Model Order")
par1.set_yticks(N)
color1 = plt.cm.viridis(0)
host.semilogy(Freq, frf_, color=color1, label="FRF")

for w_n, N in zip(nat_freqs,N):
    N= [N for i in range(len(w_n))]
    par1.plot(w_n, N,"+", markersize=5, label="Model Order")

plt.show()
#%%




































# Adjust spacings w.r.t. figsize
fig.tight_layout()

plt.savefig("pyplot_multiple_y-axis.pdf")




    plt.axvline(x = points, ls = '--' )
    
plt.show()

    #reconstruct the frf
    num = np.poly1d(a)
    denum = np.poly1d(b)
    frf_cre = num/denum
    
    s =[complex(0,1)*i for i in Freq]
    num = num(s)
    denum = denum(s)
    frf_cre = num/denum
gh =[i for i in nat_freq]

plt.semilogy(Freq,frf_)
plt.semilogy(Freq, abs(frf_cre))
