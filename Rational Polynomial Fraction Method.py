import os,pickle,zipfile
import pandas as pd
import numpy as np
import itertools

#Import and process FRF data
iters = [1]
reps = [1]
test_series = "BR_AR"

#empty list for individual sensor plot
sensor_frf_lists = {}
sensor_frf_freq_lists = {}
sensor_frf_mean = {}

for iter_ in iters:
    
    
    
    for rep in reps: 
        test_rep = rep
        test_name = test_series + '_' + str(iter_) + '_' + str(test_rep)
        zip_file_path = "Test_Data/" + test_series + "/"+ test_series \
            + "_" + str(iter_)+ '/' + test_series + "_" \
                + str(iter_) + "_" + str(rep) + '.zip'
        #directory for dumping the extraction file
        ext_dump_dir = 'Extraction_Dump'
        if not os.path.exists(ext_dump_dir):
             os.makedirs(ext_dump_dir)
        #unzip files 
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(ext_dump_dir)

        #create a list of files
        list_files = os.listdir(ext_dump_dir)
        list_files.remove('meta_data.json')
        
        #Remove FRC which is an empty list
        list_files.remove('FRC.pickle')
        
        #Data dictionary label
        data_labels = [label[:-7] for label in list_files]
        data_labels = [label.replace("-", "_" ) for label in data_labels]
        
        #create an empty list for frf and coh and their frequencies for each sensor
        if not sensor_frf_lists:
            sensor_frf_lists = {key:[] for key in list_files}
            sensor_frf_freq_lists = {key:[] for key in list_files}

        #create combine all repetition for each sensor
        for sensor in list_files:
            #open file
            infile = open(os.path.join(ext_dump_dir,sensor), 'rb')
            sensor_ = pickle.load(infile)
            infile.close()
            
            all_data = sensor_['data']
            #extract individual data
            frf= all_data[-1]
            
            #zip all data with corresponding frequencies
            frf_data, frf_freq = frf['data_y'],frf['data_x']
            
            #Append the zipped data in a list for each sensors
            sensor_frf_lists[sensor].append(frf_data)
            sensor_frf_freq_lists[sensor].append(frf_freq) 

            
    # Empty dic. for all iteration of individual sensor
    if not sensor_frf_mean:
        sensor_frf_mean = {key:[] for key in data_labels}
        sensor_frf_freq_mean = {key:[] for key in data_labels}
        
    for sensors in sensor_frf_lists:
        freqs_frf = sensor_frf_freq_lists[sensors]
        frf = sensor_frf_lists[sensors]
        
        #find average of all iterations
        frf_average = [sum(frf_list) / len(frf_list) for frf_list in zip(*sensor_frf_lists[sensors])]
        frf_freq_average = [sum(frf_freq_list) / len(frf_freq_list) for frf_freq_list in zip(*sensor_frf_freq_lists[sensors])]
        sensor = sensors[:-7]
        sensor = sensor.replace('-','_')
        #Store average for this sensor in the empty dictionary
        sensor_frf_mean[sensor] = frf_average
        sensor_frf_freq_mean[sensor] = frf_freq_average

def RFPM_alg(FRF, Freq, min_freq, max_freq, sensor, N):
    import numpy as np
    import itertools
    from scipy import signal
    """This function computes the modal parameter of a system 
    using the Rational Fraction Polynomial Method"""
    #FRF- Response Frequency as a dictionary with keys as sensor 
    #name and values as list of FRF
    #Frequency value as list 
    #Sensor input must be a the string of the sensor name
    #N - degree of freedom
    n = 2*N    #number of denominator polynomial terms
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


def Stabilization_Dia(FRF, FRF_est, Freq, Order, sensor):
    import matplotlib.pyplot as plt 
    # Create figure and subplot
    fig, host = plt.subplots(figsize=(15,10))   
    plt.grid()  
    par1 = host.twinx()
    host.set_xlabel("Frequency[Hz]")
    host.set_ylabel("FRF")
    par1.set_ylabel("Model Order")
    par1.set_yticks(order)
    color1 = plt.cm.viridis(0)
    FRF = [abs(frf) for frf in FRF]
    host.semilogy(Freq, FRF, color=color1, label="FRF")
    for i, j in zip(FRF_est, Order):
        k = [abs(frf) for frf in i]
        host.semilogy(Freq, k)
    host.legend(Order, loc="upper left")
    for w_n, N in zip(nat_freqs,Order):
        N= [N for i in range(len(w_n))]
        par1.plot(w_n, N,"+", markersize=5, label='gaga')
    plt.title(sensor)
    plt.show()
    return fig

sensor = 'EXH'
N = [i for i in range(2,10)]
nat_freqs =[]
damp_ratio =[]
order = []
frf_est = []
for i in N:
    nat_freq_, dp, order_, FRF, Freq, FRF_est_ = RFPM_alg(sensor_frf_mean, sensor_frf_freq_mean, 20, 50, sensor, i)
    nat_freqs.append(nat_freq_)
    order.append(order_)
    damp_ratio.append(dp)
    frf_est.append(FRF_est_)
    
plot = Stabilization_Dia(FRF,frf_est, Freq, order, sensor)


