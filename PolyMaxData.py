def PolyMaxDataPrep(data):
    """This function is pecific to the data sent available, it coverts the two dimensional 
    space data (SIMO) which is a dictionary with keys equal to the sensor (output)
    and vales equal to a list of each input across all frequencies to a three dimensional
    matix space
    e.g [[[FRF(w) for input 1 for ouput 1 ], [FRF(w) for input 2 for ouput 1 ], ...],
         [[FRF(w) for input 1 for ouput 2 ], [FRF(w) for input 2 for ouput 2 ], ...]], ...]"""
    import pandas as pd
    import numpy as np
    #Converting file to data frame and array
    data = pd.DataFrame(data)
    data = data.to_numpy()
    #Converting data from 2D [Nf*l] to 3D [m*NF*l]
    data = np.expand_dims(data, axis= 0)
    return data
# =============================================================================
#     DATA = []
#     for i in data.values():
#         dump_dim = []
#         dump_dim.append(i)
#         DATA.append(dump_dim)
#     return DATA
# =============================================================================
