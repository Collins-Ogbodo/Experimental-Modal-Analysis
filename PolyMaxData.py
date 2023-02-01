def PolyMaxDataPrep(data):
    """This function is pecific to the data sent available, it coverts the two dimensional 
    space data (SIMO) which is a dictionary with keys equal to the sensor (output)
    and vales equal to a list of each input across all frequencies to a three dimensional
    matix space
    e.g [[[FRF(w) for input 1 for ouput 1 ], [FRF(w) for input 2 for ouput 1 ], ...],
         [[FRF(w) for input 1 for ouput 2 ], [FRF(w) for input 2 for ouput 2 ], ...]], ...]"""
    DATA = []
    for i in data.values():
        dump_dim = []
        dump_dim.append(i)
        DATA.append(dump_dim)
    return DATA