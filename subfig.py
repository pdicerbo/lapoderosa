import numpy as np
import os
import string
from bisect import bisect_left # for BilinearInterpolation
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
import colorsys
#from matplotlib.mlab import griddata

matrix_Logdelta_LogT_H2       = 'matrix_modif_Logdelta_LogT_H2.dat'
matrix_Logdelta_LogT_H2_tcool = 'matrix_modif_Logdelta_LogT_tcool.dat'
path_out                      = '/scratch2/dicerbo/plot_test/prova/'
path_plot                     = '/scratch2/dicerbo/plot_test/newexit/'

# global arrays: Temperature, H2OverDensity, H2Fraction, tcool to load UM's tables
#                T in K, tcool in Gyr
T          = None          # dimension 1x50
Dens       = None          # dimension 1x50
FH2        = None          # dimension 50x50
t_cool     = None          # dimension 50x50


def main():
    global path_out
    '''
    directory = '/home/dicerbo/output/scratch2/plot_test/prova/prova1'
    if not os.path.exists(directory):
        os.makedirs(directory)
        print '%s created successfully'%(directory)
    '''
    dirs = os.listdir(path_out)
    for d in dirs:
        if string.count(d, 'jpg') == 0 and string.count(d, 't50') == 1:
            print '\n\tStart working on '+ d
            #adjust(path_out, d)
            plot_def(d)
            print '\n\tEnd working on ' + d

    print '\n\tFinally End\n'

def adjust(path, directory):

    print '\n\tWithin adjust function\n'

    global matrix_Logdelta_LogT_H2
    LoadMatrix(filename=matrix_Logdelta_LogT_H2)
    global T ; global Dens ; global FH2

    files = os.listdir(path+directory)
    for name in files:
        if string.count(name, 'time') != 0:
            #print '\tWorking on '+name
            matrix = np.loadtxt(path+directory+'/'+name, comments = '#')
            if matrix.size <= 6:
                print '\n\t'+name+' has '+str(matrix.size / 3)+' lines\n'
                if matrix.size == 3:
                    mt = np.zeros((3, 3), dtype = float)
                    mt[0:] = matrix[0:]
                    fp = open(path+directory+'/'+name, 'w')
                    np.savetxt(fp, mt, fmt='%g', delimiter='\t', newline='\n')
                    fp.flush(); fp.close()
                    print '\t-> Corrected\n'
                elif matrix.size == 6:
                    mt = np.zeros((3, 3), dtype = float)
                    mt[0] = matrix[0]
                    mt[1] = matrix[0]
                    mt[2] = matrix[1]
                    fp = open(path+directory+'/'+name, 'w')
                    np.savetxt(fp, mt, fmt='%g', delimiter='\t', newline='\n')
                    fp.flush();fp.close()
                    print '\t-> Corrected\n'
                else:
                    print '\n\tFile'+name+' have some problem.....\n'

            else:
                print '\tFile '+name+' have enough data'
        else:
            print "\n\tFile " + name + " is for Blitz&Rosolowsky's plot -> Continue\n"
            continue


def plot_def(directory):
    print '\n\tWithin plot function\n'
    global path_read; global path_out
    #Load tcool matrix
    LoadMatrix(filename=matrix_Logdelta_LogT_H2_tcool)
    global T ; global Dens ; global FH2; global t_cool

    tcool = t_cool
    tcool[tcool > 0.] = np.log10(tcool[tcool > 0.])
    v_min = -5
    v_max = 7.
    tcool[tcool == 0.] = v_min
    tcool[tcool > v_max] = v_max
    tcool[tcool <= v_min] = v_max
    '''
    H2 = FH2
    H2[H2 > 0.] = np.log10(H2[H2 > 0.])
    v_min = -6
    v_max = -2.
    H2[H2 == 0.] = v_min
    H2[H2 > v_max] = v_max
    H2[H2 < v_min] = v_min
    '''
    numlev = 15
    dmag0 = (v_max - v_min) / float(numlev)
    levels0 = np.arange(numlev) * dmag0 + v_min
    
    #path's plot
    files = os.listdir(path_out+directory)
    files.sort()
    fls = files[:]
    press = np.zeros(len(files), dtype = float)
    j = 0
    for name in files:
        if string.count(name, 'time') != 0:
            fls[j] = directory+'/'+name
            press[j] = float(name[(len(name)-7):-4])
            j += 1
        else:
            br = path_out + directory + '/' + name
            print "\n\tFile " + name + " is for Blitz&Rosolowsky's plot -> Continue\n"

    if j == len(files):
        filedef = fls[:]
        pdef = press[:]
    else:
        filedef = fls[:(j-len(files))]
        pdef = press[:(j-len(files))]
    pmax = pdef.max()
    pmin = pdef.min()

    h = np.zeros(pdef.size, dtype = float)
    ind = 0
    for p in pdef:
        h[ind] = ((p-pmin) / (pmax-pmin))*250.
        ind += 1
    cdef = [colorsys.hsv_to_rgb(x/360., 1., 1.) for x in h]
    
    #plots
    fig = plt.figure(figsize=(18,16))
    figura = fig.add_subplot(2, 1, 1, adjustable='box', aspect = 1.1)
    plt.title('Paths in Phase Diagram', fontsize = 28)
    #figura = plt.contourf(Dens,T,H2,levels0,extend='both', cmap = cm.hot)
    figura = plt.contourf(Dens,T,tcool,levels0,extend='both', cmap = cm.hot_r)
    ax1 = plt.gca()
    ax1.set_xlim([Dens.min(), Dens.max()])
    ax1.set_ylim([1., 5.])
    cbar = plt.colorbar(figura,format='%3.1f', shrink=0.8)
    cbar.set_ticks(np.linspace(v_min,v_max,num=levels0.size,endpoint=True))
    #cbar.set_label('H$_{2}$ fraction',fontsize=20)
    cbar.set_label('t_cool',fontsize=20)
    print "\n\tUmberto's matrix plotted\n"
        
    k = 0
    for name in filedef:
        print '\tPlotting ' + name[(len(directory)+1):] + ' file'
        #figura = plt.plotfile(path_out+name, delimiter = '\t', cols=(1, 2), comments='#', color = cdef[k], marker='.', 
                                         #mfc = cdef[k], mec = cdef[k], label = 'Log10P = '+str(pdef[k]), newfig=False)
        data = np.loadtxt(path_out+name, comments = '#'); data = data.T
        rho = data[1, :]; tmp = data[2, :]
        plt.plot(rho, tmp, color = cdef[k], marker='.', mfc = cdef[k], mec = cdef[k], label = 'Log10P = '+str(pdef[k]))
        k += 1
    lgd = plt.legend(bbox_to_anchor=(1.75, 0.4), loc=5, borderaxespad=1.)
    
    plt.xlabel('log10 $\\rho$',fontsize=20) ; plt.ylabel('Log10 T[k]',fontsize=20)

    #Blitz&Rosolowsky plot
    figura2 = fig.add_subplot(2, 1, 2, adjustable='box', aspect = 1.3)
    plt.title('Blitz & Rosolowsky', fontsize = 28)
    ax2 = plt.gca()
    newm = np.loadtxt(br, comments = '#'); newm = newm.T
    press = newm[0, :]
    br_ro = newm[3, :]
    fh2   = newm[4, :]
    ax2.set_xlim([3., 6.])
    ax2.set_ylim([0., 1.02])
    plt.plot(press, br_ro, 'k-')
    plt.plot(press, fh2, 'b-')
    #scale figure2
    scale = figura2.get_position().bounds
    newpos = [scale[0]*3./4. + 0.2, scale[1]*3./4., scale[2]*3./4., scale[3]*3./4.]
    figura2.set_position(newpos)

    newname = path_plot + 'path_' + directory + '.jpg'
    plt.savefig(newname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    #plt.savefig(newname)
    plt.close('all')
    print '\n\t'+newname[len(path_plot):]+' done\n'

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
