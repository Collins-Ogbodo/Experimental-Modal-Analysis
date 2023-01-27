#Import the needed packages
import os,pickle,zipfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
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
#%%       
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
    #M[:,-1] = 1.0
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
    return w_n, M, CM, R_sensor, M_sensor
    
#%%
#Convert average frf for all sensor to dataframe
data_frf = pd.DataFrame(sensor_frf_mean) 
FRF = data_frf
FRF['mean_frf'] = FRF.mean(axis=1)
#converting dataframe to array of dataframe 
data_frf = data_frf.to_numpy()
frf = data_frf
FRF = FRF.to_numpy()
FRF = FRF[:,-1]
#Extracting frequency for a sensor since frequency values for all sensors are equal
Freq = sensor_frf_freq_mean['EXH']

#Frequency range under consideration
#Freq = Freq[0:800]
#frf =  frf[0:800]
#FRF = FRF[0:800]
N = [i for i in range(2,10)]
nat_freqs =[]

order = []
for i in N:
    W_n, M, MC, R, M = PolyMAX(frf, Freq, i)
    #nat_freq_, dam_ratio_, b_, poles_, order_ = G_RFPM_alg(frf, Freq, i)
    nat_freqs.append(W_n)
    order.append(i)

#Extract the amplitufe of the frf
frf_ = [abs(i) for i in FRF]

# Create figure and subplot
fig, host = plt.subplots(figsize=(15,10))   
plt.grid()  
par1 = host.twinx()
host.set_xlabel("Frequency[Hz]")
host.set_ylabel("FRF")
par1.set_ylabel("Model Order")
par1.set_yticks(order)
color1 = plt.cm.viridis(0)
host.semilogy(Freq, frf_, color=color1, label="FRF")

for w_n, N in zip(nat_freqs,order):
    N= [N for i in range(len(w_n))]
    par1.plot(w_n, N,"+", markersize=5, label="Model Order")
plt.title('Global-RFPM')
plt.show()


#%%
W_n, M, MC, R, M = PolyMAX(frf, Freq, 20)

#%%






