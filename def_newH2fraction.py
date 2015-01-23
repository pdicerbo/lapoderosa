import numpy as np
import os
from bisect import bisect_left # for BilinearInterpolation
import matplotlib.pyplot as plt

matrix_Logdelta_LogT_H2       = 'matrix_modif_Logdelta_LogT_H2.dat'
matrix_Logdelta_LogT_H2_tcool = 'matrix_modif_Logdelta_LogT_tcool.dat'
path_out                      = '/scratch2/dicerbo/plot_path/sim_def/'
# global variables
redshift          = 19.0
Hubble            = 0.72
h2                = Hubble * Hubble
Omega0m           = 0.3
Omega0l           = 0.7
Omega0r           = 0.0
BOLTZMANN         = 1.3806e-16  # erg/K
PROTONMASS        = 1.6726e-24  # g
GAMMA_MINUS1      = 5./3. - 1.
HYDROGEN_MASSFRAC = 0.76
mu_h              = 4./ (5. * HYDROGEN_MASSFRAC + 3.)  #  molecular weight of the hot phase
mu_c              = 4./ (3. * HYDROGEN_MASSFRAC + 1.)  # molecular weight of the cold phase
FracC             = 0.9
NPCLASS           = 300
rho_cr            = 1.9e-29 * ((1-Omega0m-Omega0l)*pow((1+redshift),2) + Omega0m*pow((1+redshift),3) + Omega0r*pow((1+redshift),4) + Omega0l) * h2


UnitMass_in_g            = 1.989e43
UnitLength_in_cm         = 3.085678e21
UnitVelocity_in_cm_per_s = 1.0e5
UnitTime_in_s            = UnitLength_in_cm / UnitVelocity_in_cm_per_s
UnitEnergy_in_cgs        = UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2)
UnitDensity_in_cgs       = UnitMass_in_g / pow(UnitLength_in_cm,3)

SE_TO_T                  = GAMMA_MINUS1 * UnitEnergy_in_cgs * mu_h * PROTONMASS / BOLTZMANN / UnitMass_in_g
T_TO_SE                  = 1.0/SE_TO_T
Th2                      = 10.0
FTHR                     = 1.e-6

# global arrays: Temperature, H2OverDensity, H2Fraction, tcool to load UM's tables
#                T in K, tcool in Gyr
T          = None          # dimension 1x50
Dens       = None          # dimension 1x50
FH2        = None          # dimension 50x50
t_cool     = None          # dimension 50x50

def main():
    # functions to call
    options = {1:'OneT',2:'OneDelta',3:'Fit_T_Delta',4:'MolecularProfileI',5:'MolecularProfileTc',6:'plot'}
    num=1
    while num in range(1,7):
        # select what you want to do...
        print '\n\t    Select                    Action'
        print   '\t    1      H2 fraction as a function of the overdensity at a fix T'
        print   '\t    2      H2 fraction as a function of the temperature at a fixed overdensity'
        print   '\t    3      H2 fraction given one temperature and one overdensity (fitted)'
        print   '\t    4      H2 molecular profile as a function of the pressure, given an SPH density (istantaneous)'
        print   '\t    5      H2 molecular profile as a function of the pressure, given an SPH density (converging Tc)'
        print   '\t    6      make the plots'
        print   '\t    whatever else             Exit'
        #num = int(raw_input("\n\t Select Action :> "))
        num = 5
        if num in options.keys():
            eval(options[num])()
            num = 10

    print '\n\t END OF GAME!!!\n'


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



def plot():
    """
    This function makes the following plots:
    - T vs Density and as color the H2Fraction
    - T vs Density and as color the cooling time
    """

    # first plot
    global matrix_Logdelta_LogT_H2
    LoadMatrix(filename=matrix_Logdelta_LogT_H2)
    global T ; global Dens ; global FH2
    H2 = FH2
    H2[H2 > 0.] = np.log10(H2[H2 > 0.])
    vmin = -6
    vmax = -2.
    H2[H2 == 0.] = vmin
    H2[H2 > vmax] = vmax
    H2[H2 < vmin] = vmin
    nlev = 15
    dmag = (vmax - vmin) / float(nlev)
    levels = np.arange(nlev) * dmag + vmin
    plt.figure()
    figura = plt.contourf(Dens,T,H2,levels,extend='both')
    figura0 = plt.contour(Dens,T,H2,levels,colors = ('k',),linewidths = (0.3,))
    plt.title('H$_{2}$ fraction')
    plt.xlabel('log10 $\delta$',fontsize=20) ; plt.ylabel('Log10 T[k]',fontsize=20)
    cbar = plt.colorbar(figura,format='%3.1f')
    cbar.set_ticks(np.linspace(vmin,vmax,num=levels.size,endpoint=True))
    cbar.set_label('H$_{2}$ fraction',fontsize=20)
    plt.savefig('fraction_H2.jpg')
    plt.close('all')
    print '\n\t fraction_H2.jpg done\n'

    # second plot
    global matrix_Logdelta_LogT_H2_tcool
    LoadMatrix(filename=matrix_Logdelta_LogT_H2_tcool)
    global t_cool
    cool = t_cool
    cool[cool > 0.] = np.log10(cool[cool > 0.])
    vmin = -5
    vmax = 7.
    cool[cool == 0.] = vmin
    cool[cool > vmax] = vmax
    cool[cool < vmin] = vmin
    nlev = 15
    dmag = (vmax - vmin) / float(nlev)
    levels = np.arange(nlev) * dmag + vmin
    plt.figure()
    figura = plt.contourf(Dens,T,cool,levels,extend='both')
    figura0 = plt.contour(Dens,T,cool,levels,colors = ('k',),linewidths = (0.3,))
    plt.title('H$_{2}$ fraction')
    plt.xlabel('log10 $\delta$',fontsize=20) ; plt.ylabel('Log10 T[k]',fontsize=20)
    cbar = plt.colorbar(figura,format='%3.1f')
    cbar.set_ticks(np.linspace(vmin,vmax,num=levels.size,endpoint=True))
    cbar.set_label('tcool   [Gyr]',fontsize=20)
    plt.savefig('fraction_H2_tcool.jpg')
    plt.close('all')
    print '\n\t fraction_H2_tcool.jpg done\n'
    
    
def OneT():
    """
    H2 fraction as a function of the overdensity at a fix T
    """

    print '\n\t Within OneT function\n\n'

    global T ; global Dens ; global FH2
    if T==None or Dens==None or FH2==None:
        # load the matrix
        LoadMatrix(filename=matrix_Logdelta_LogT_H2)
        
    for i in range(T.size):
        print 'index:{0:2d}    logT:{1:f}    T:{2:e}'.format(i,T[i],10.0**T[i])

    choice = int(raw_input('\n\t Enter T value (integer)>'))
    # select only FH2 at this Temperature
    H2 = FH2[choice,:]
    filename = 'oneT'+str(choice)+'.dat'
    header = '# logT:{0:f}   T:{1:e}'.format(T[choice],10.0**T[choice])
    header += '\n# 1] Overdensity'
    header += '\n# 2] H2 fraction\n#\n'
    fp = open(filename,'w')
    fp.write(header)
    np.savetxt(fp,np.array([Dens,H2]).T,fmt='%e   %e',delimiter='   ',newline='\n')
    fp.flush() ; fp.close()
    print '\n\t Output in file '+filename+'\n\n'


def OneDelta():
    """
    H2 fraction as a function of the temperature at a fixed overdensity
    """
    
    print '\n\t Within OneDelta function\n\n'

    global T ; global Dens ; global FH2
    if T==None or Dens==None or FH2==None:
        # load the matrix
        LoadMatrix(filename=matrix_Logdelta_LogT_H2)

    for i in range(Dens.size):
        print 'index:{0:2d}    logDens:{1:f}    Dens:{2:e}'.format(i,Dens[i],10.0**Dens[i])

    choice = int(raw_input('\n\t Enter Density value (integer)>'))
    # select only FH2 at this OverDensity
    H2 = FH2[:,choice]
    filename = 'oneDelta'+str(choice)+'.dat'
    header = '# logDens:{0:f}   Dens:{1:e}'.format(Dens[choice],10.0**Dens[choice])
    header += '\n# 1] Temperature'
    header += '\n# 2] H2 fraction\n#\n'
    fp = open(filename,'w')
    fp.write(header)
    np.savetxt(fp,np.array([T,H2]).T,fmt='%e   %e',delimiter='   ',newline='\n')
    fp.flush() ; fp.close()
    print '\n\t Output in file '+filename+'\n\n'
    
def Fit_T_Delta():
    """
    This function returns the H2 fraction given one temperature and one overdensity (fitted)
    """
    
    print '\n\t Within Fit_T_Delta function\n\n'
    
    global T ; global Dens ; global FH2
    if T==None or Dens==None or FH2==None:
        # load the matrix
        LoadMatrix(filename=matrix_Logdelta_LogT_H2)

    Tmin = T.min() ; Tmax = T.max()
    Densmin = Dens.min() ; Densmax = Dens.max()

    print '\t T  available range: Tmin  {0:e} and Tmax  {1:e}'.format(10.**Tmin,10.**Tmax)
    print '\t Density available range: Densmin {0:e} and Densmax {1:e}'.format(10.**Densmin,10.**Densmax)

    while True:
        Temp, Density = raw_input('\n\tEnter Temperature and Density (LINEAR)>').split()
        #Temp = np.log10(float(Temp)) ; OverDens = np.log10(float(OverDens))
        if Temp>=Tmin and Temp<=Tmax and Density>=Densmin and Density<=Densmax: break

    H2FractionInterpolated = float(BilinearInterpolation(Dens, T, FH2, Density, Temp))
    print '\n\t Temp: {0:e}, Density: {1:e}, InterpFracH2: {2:e}\n'.format(10.**Temp,10.**Density,H2FractionInterpolated)
    
def FitTDelta(Temp, Density, matrix):
    """
    This function returns the H2 fraction given one temperature and one overdensity (fitted)
    """
    global T ; global Dens 
    H2FractionInterpolated = float(BilinearInterpolation(Dens, T, matrix, Density, Temp))
    #print '\n\t Temp: {0:e}, OverDensity: {1:e}, InterpFracH2: {2:e}\n'.format(10.**Temp,10.**OverDens,H2FractionInterpolated)
    return H2FractionInterpolated

def BilinearInterpolation(x_values, y_values, values, x, y):
    """
    This function performs the Bilinear Interpolation

    x_values is array Overdensity given by Umberto's matrix
    y_values is array Temperature given by Umberto's matrix
    values   is grid  H2 fraction given by Umberto's matrix
    x        is the enter overdensity
    y        is the enter temperature

    N.B. This function does not check boundary conditions yet

    """
    x = np.log10(x)
    y = np.log10(y)
    j = bisect_left(x_values, x) - 1
    i = bisect_left(y_values, y) - 1
    #condizioni di uscita PROVVISORIE
    if ((x_values.max() - x) <= 0. or (x_values.min() - x) >= 0.) and ((y_values.max() - y) <= 0 or (y_values.min() - y) >= 0): 
        return 0.
    elif ((x_values.max() - x) >= 0. and (x_values.min() - x) <= 0.) and ((y_values.max() - y) >= 0 and (y_values.min() - y) <= 0):
        #tipical case
        x1 = x_values[j:j + 1]
        x2 = x_values[j + 1:j + 2]
        y1 = y_values[i:i + 1] 
        y2 = y_values[i + 1:i + 2]
        f1, f3 = values[i][j:j + 2]
        f2, f4 = values[i + 1][j:j+2]
    
        dx = x - x1; dy = y - y1
        lx = x2 - x1; ly = y2 - y1

    elif x_values.max() < x:
        x1 = x_values[j - 1:j]
        x2 = x_values[j:j + 1]
        y1 = y_values[i:i + 1]
        y2 = y_values[i + 1:i + 2]
        f1 = values[i][j]; f3 = 0.
        f2 = values[i + 1][j]; f4 = 0.
        
        dx = 0.; dy = y - y1
        lx = x2 - x1; ly = y2 - y1

    elif x_values.min() > x:
        x1 = x_values[j:j + 1]
        x2 = x_values[j + 1:j + 2]
        y1 = y_values[i:i + 1]
        y2 = y_values[i + 1:i + 2]
        f1 = 0.; f3 = values[i][j]
        f2 = 0.; f4 = values[i + 1][j]

        dx = 0.; dy = y - y1
        lx = x2 - x1; ly = y2 - y1

    elif y_values.min() > y:
        x1 = x_values[j:j + 1]
        x2 = x_values[j + 1:j + 2]
        y1 = y_values[i:i + 1]
        y2 = y_values[i + 1:i + 2]
        f1 = 0.; f3 = 0.
        f2 = values[i][j]; f4 = values[i][j + 1]

        dx = x - x1; dy = 0.
        lx = x2 - x1; ly = y2 - y1

    elif y_values.max() < y:
        x1 = x_values[j:j + 1]
        x2 = x_values[j + 1:j + 2]
        y1 = y_values[i - 1:i]
        y2 = y_values[i:i + 1]
        f1 = values[i][j]; f3 = values[i][j + 1]
        f2 = 0.; f4 = 0.;

        dx = x - x1; dy = 0.
        lx = x2 - x1; ly = y2 - y1

    f_fit = f1 + (f3 - f1)*dx/lx + (f2 - f1)*dy/ly + (f1 + f4 - f2 - f3)*dx*dy/(lx*ly)

    return f_fit

def molecular_fraction(pressure):
    """
    This function takes the pressure as input
    and returns the molecular fraction following
    Blitz & Rosolowski (2006) relation
    N.B. pressure must be physical (no h100)
    """

    PRESFMOL = 20000.

    return  1. / (1. + (PRESFMOL/pressure))

def MolecularProfileI():
    print '\t\n Within MolecularProfileI function\n'
    
    global T ; global Dens ; global FH2 ;
    global BOLTZMANN; global GAMMA_MINUS1; global mu_c; #devo richiamare anche le costanti come variabili globali?
    if T==None or Dens==None or FH2==None:
        # load the matrix                                                                                                                                                               
        LoadMatrix(filename=matrix_Logdelta_LogT_H2)
    
    T_c = 2500.;
    E_c = BOLTZMANN * T_c/GAMMA_MINUS1/mu_c/PROTONMASS;
    
    Pmass=2.78e-4; #GA1 initial mass
    Density=0.05; #in the middle of SF MUPPI cloud in phase diagram
    
    #Pmin & Pmax definition
    pmin = 1.e2
    pmax = 1.e6
    
    #logaritmic interval
    dp = (np.log10(pmax) - np.log10(pmin))/NPCLASS
    
    #define list for cycle
    lst = np.arange(np.log10(pmin), np.log10(pmax), dp)

    filename = 'fh2-vs-press.Tc'+str(T_c)+'.dat'
    header = '#Press\tRho_c\tFcoll\tFH2'
    header += '\n#\n'
    
    fp = open(filename,'w')
    fp.write(header)
    #np.savetxt(fp,np.array([T,H2]).T,fmt='%e   %e',delimiter='   ',newline='\n')
    
    for p in lst:
        press = 10.**p
        Rho_c = press/(GAMMA_MINUS1*E_c)*BOLTZMANN #in cgs
        #OverDensity = Rho_c/rho_cr
        Fcoll = molecular_fraction(press)
        frac = float(BilinearInterpolation(Dens, T, FH2, Rho_c, T_c))
        #line to write
        #line = str(press)+'\t'+str(Rho_c)+'\t'+str(Fcoll)+'\t'+str(frac)+'\n'
        line = '%g\t%e\t%e\t%e\n'%(p, Rho_c, Fcoll, frac)
        fp.write(line)

    fp.flush() ; fp.close()
    print '\n\t Output in file '+filename+'\n\n'



FracH = 0.
t_atomic = 0.
dens_a = 0.
mmol = 0.
index = 0
filectrl = 0

def MolecularProfileTc():
    print '\n\t Within MolecularProfileTc function\n'
    
    global T; global Dens; global FH2; global t_cool
    global mmol; global t_atomic; global dens_a; global index
    if T==None or Dens==None or FH2==None:
        # load the matrix                                                                                                                                                               
        LoadMatrix(filename=matrix_Logdelta_LogT_H2)
    if t_cool==None:
        LoadMatrix(filename=matrix_Logdelta_LogT_H2_tcool)

    fh2 = FH2
    tcool = t_cool*1.e9    #tcool in yr
    max_tcool = tcool.max()
    tcool[tcool <= 1.e2] = max_tcool

    Pmass = 2.78e-4; #GA1 initial mass
    Density = 0.05; #in the middle of SF MUPPI cloud in phase diagram
    T_a = 900.
    if os.path.exists(path_out+'t'+str(T_a)):
        print 'path %s exist'%('t'+str(T_a))
    else:
        to_make = path_out+'t%s'%(str(T_a))
        os.makedirs(to_make[:-2])
        print '%s created successfully'%('t'+str(T_a))
    newpath = to_make[:-2]+'/'
    pmin = 1.e3
    pmax = 1.e6
    dp = (np.log10(pmax) - np.log10(pmin))/NPCLASS
    lst = np.arange(np.log10(pmin), np.log10(pmax), dp)

    #time = float(raw_input("\n\t Enter total time [Gyr] :> "))
    #time = time*10**9
    time = 1.e9    #total time of simulation
    tstep = 4.*1.e3   #timestep (smaller then minimum of tcool matrix)
    
    filename = newpath+'press-vs-fracH2.T_atom'+str(T_a)+'.dat'
    header = '#Press\tRho_atom\tT_atom\tFcoll\tFH2\tM_H2'
    header += '\n#\n'

    fp = open(filename, 'w')
    fp.write(header)
    #print '%s ready to write'%(filename[len(newpath):])
    pstr = 0
    for p in lst:
        press = 10.**p   #P/k_b
        if press*mu_c*PROTONMASS/T_a < 10.**Dens.min() or press*mu_c*PROTONMASS/T_a > 10.**Dens.max():
            if pstr == 0:
                pstart = p; tstart = T_a; pstr = 10.
            Fcoll = molecular_fraction(press)
            frac = 0.
            Rho_a = press*mu_c*PROTONMASS/T_a
            line = '%g\t%e\t%g\t%e\t%e\t%e\n'%(p, dens_a, t_atomic, Fcoll, frac, mmol)
            fp.write(line)
            fp.flush()
            index += 1
        else:
            Fcoll = molecular_fraction(press)
            frac = time_int(press, Pmass, Density, T_a, tcool, fh2, time, tstep, newpath)
            Rho_a = press*mu_c*PROTONMASS/T_a
            line = '%g\t%e\t%g\t%e\t%e\t%e\n'%(p, dens_a, t_atomic, Fcoll, frac, mmol)
            fp.write(line)
            fp.flush()
            index += 1
    if pstr > 1:
        line = '#(Log10)Pressure %g at temperature %g is out of bound\n'%(pstart, tstart)
        fp.write(line)
    fp.write('# Done!')
    fp.flush()
    fp.close()

    print '\n\t Output in file '+filename+'\n\n'


def time_int(press, Pmass, Density, T_a, tcool, fh2, time, tstep, path):
    #global t_cool
    global Dens
    global t_atomic, dens_a
    global mmol
    global index
    global filectrl
    Temp = T_a
    mass = Pmass * UnitMass_in_g/Hubble
    ma = mass*FracC

    #tcool = t_cool*1.e9
    #tcoolmax = tcool.max()

    if index == 10:
        fname = path+'time_evolution_log10P'+str(np.log10(press))+'.dat'
        header = '# time\tLog10Rho_a\tLog10_T'
        header += '\n#\n'
        wr = open(fname, 'w')
        wr.write(header)

    Rho_a = press*mu_c*PROTONMASS/Temp
    mh2 = 0.
    dmin = Dens.min(); dmax = Dens.max()
    dmin = 10.**dmin; dmax = 10.**dmax
    t = 0.
    ctrl = 2
    while t < time:
        if ma <= 0. or press*mu_c*PROTONMASS/Temp < dmin or press*mu_c*PROTONMASS/Temp > dmax:
            if ma <= 0.:
                if index == 10:
                    if filectrl == 0:
                        line = '%e\t%e\t%e\n'%(0,np.log10(Rho_a), np.log10(Temp))
                        line += '%e\t%e\t%e\n'%(8e3,np.log10(Rho_a), np.log10(Temp))
                        line += '%e\t%e\t%e\n'%(16e3,np.log10(Rho_a), np.log10(Temp))
                        wr.write(line);
                    elif filectrl == 1:
                        line = '%e\t%e\t%e\n'%(8e3,np.log10(Rho_a), np.log10(Temp))
                        line += '%e\t%e\t%e\n'%(16e3,np.log10(Rho_a), np.log10(Temp))
                        wr.write(line);
                    elif filectrl == 2:
                        line = '%e\t%e\t%e\n'%(16e3,np.log10(Rho_a), np.log10(Temp))
                        wr.write(line);                        
                    wr.flush(); wr.close()
                    t_atomic = Temp
                    dens_a = Rho_a
                    mmol = mass*FracC
                    index = 0
                filectrl = 0
                return 1.
            elif mh2/(mh2 + ma) > 1.e-2 and np.log10(Temp) < 3.:
                if index == 10:
                    if filectrl == 0:
                        line = '%e\t%e\t%e\n'%(0, np.log10(Rho_a), np.log10(Temp))
                        line += '%e\t%e\t%e\n'%(8e3, np.log10(Rho_a), np.log10(Temp))
                        line += '%e\t%e\t%e\n'%(16e3,np.log10(Rho_a), np.log10(Temp))
                        wr.write(line);
                    elif filectrl == 1:
                        line = '%e\t%e\t%e\n'%(8e3, np.log10(Rho_a), np.log10(Temp))
                        line += '%e\t%e\t%e\n'%(16e3, np.log10(Rho_a), np.log10(Temp))
                        wr.write(line);
                    elif filectrl == 2:
                        line = '%e\t%e\t%e\n'%(16e3, np.log10(Rho_a), np.log10(Temp))
                        wr.write(line);                        
                    wr.flush(); wr.close()
                    t_atomic = Temp
                    dens_a = Rho_a
                    mmol = mh2
                    index = 0
                filectrl = 0
                return 1.
            else:
                filectrl = 0
                return mh2/(mh2 + ma)
        else: 
            frac_h2 = FitTDelta(Temp, Rho_a, fh2)
            if frac_h2 < 1.e-7:
                frac_h2 = 0.
            tcooling = FitTDelta(Temp, Rho_a, tcool)
            if tcooling <= 1.e2:
                tcooling = 1.e9    #from bilinear_interpolation can return 0.
            tfact = tfactor(tcooling, tstep)   #return tstep/tcooling
            mass_f = 2.*frac_h2/(1.-frac_h2)    #conversion factor between number density fraction and mass fraction
            mh2 += ma*mass_f*tfact
            ma -= ma*mass_f*tfact
            tmp = Temp
            if index == 10 and ctrl == 2:
                line = '%e\t%e\t%e\n'%(t,np.log10(Rho_a), np.log10(Temp))
                wr.write(line)
                filectrl += 1
                ctrl = 0
            Temp = Temperature(Temp, tcooling, tstep)   #return Temp*(1 - tstep/tcooling)
            Rho_a = Rho_a*tmp/Temp
            t += tstep
            ctrl += 1
    
    t_atomic = Temp
    dens_a = Rho_a
    mmol = mh2
    if index == 10:
        if filectrl == 0:
            line = '%e\t%e\t%e\n'%(0,np.log10(Rho_a), np.log10(Temp))
            line += '%e\t%e\t%e\n'%(8e3,np.log10(Rho_a), np.log10(Temp))
            line += '%e\t%e\t%e\n'%(16e3,np.log10(Rho_a), np.log10(Temp))
            wr.write(line);
        elif filectrl == 1:
            line = '%e\t%e\t%e\n'%(8e3,np.log10(Rho_a), np.log10(Temp))
            line += '%e\t%e\t%e\n'%(16e3,np.log10(Rho_a), np.log10(Temp))
            wr.write(line);
        elif filectrl == 2:
            line = '%e\t%e\t%e\n'%(16e3,np.log10(Rho_a), np.log10(Temp))
            wr.write(line);                        
        wr.flush()
        wr.close()
        index = 0
    filectrl = 0
    return mh2/(mh2 + ma)
    
def tfactor(tcooling, tstep):
    if tstep/tcooling <= 1.:
        return tstep/tcooling
    else:
        return 1.

def Temperature(Temp, tcooling, tstep):
    f = tfactor(tcooling, tstep)
    tmp = Temp*(1 - f)
    if tmp >= 10.:
        return tmp
    else:
        return 10.
main()
