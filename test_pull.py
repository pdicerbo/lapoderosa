import numpy as np
import os
import string
from bisect import bisect_left # for BilinearInterpolation
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
import colorsys

def LoadMatrix(filename=False):
    """
    This function loads one Umberto's file,
    returns the matrix and the corresponding edges

    """

    global matrix_Logdelta_LogT_H2
    global matrix_Logdelta_LogT_H2_tcool

    if filename==False:
        raise IOError('\n\t filename MUST be provided\n')
    
    # store the path of this module
    # locate = inspect.getfile(LoadMatrix)
    # dir_file = locate.replace('H2fraction.py','')
    # filename = dir_file+filename
    if not os.path.isfile(filename):
        raise IOError('\n\t filename ',filename,' NOT found\n')

    # load file
    matrix = np.loadtxt(filename,comments='#')

    # OverDensity edges
    global Dens ; global T ; global FH2 ; global t_cool
    Dens = matrix[0,:]
    # Temperature edges
    T = matrix[1,:]

    if filename == matrix_Logdelta_LogT_H2:
        FH2 = matrix[2:,:]
    elif filename == matrix_Logdelta_LogT_H2_tcool:
        t_cool = matrix[2:,:]
    else:
        raise IOError('\n\t It seems that ',filename,' does not exist\n')


main()
