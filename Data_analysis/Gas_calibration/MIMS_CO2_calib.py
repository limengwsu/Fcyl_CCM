#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 14:33:59 2020

@author: LIMeng
"""

#import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages#save multi fig to pdf
from scipy.optimize import curve_fit
from scipy import stats
from scipy import interpolate
from scipy import ndimage as ndi

def show_plot():
    plt.legend()
    plt.show()
    plt.close()
def I_to_P(Igas, k0, k1, k2, k3, k4):#k4 for all co2
    #Igas in 1e-14 A
    #pressure pH2O, P in 1e-14 mbar, or 1e-12Pa, pPa
    #k values in 0.1A/mbar, or mA/Pa
    k_IP = np.array([[k0, k1, k2, k3, k4, k4]], dtype=float)
    Pgas = Igas / k_IP#Igas should exclude H2O
    #Igas should be a (n,x) shpaed array, recording the current
    #k_IP = I/P, also a (1,x) shaped array for calculations
    P = np.sum(Pgas, axis= 1)
    #pH2O includes H2O and other gas, treat as constant
    return P
def pgas(Igas, paras):
    #paras is the popt result from curve_fit, extract k1 to k4 for calculation
    k0, k1, k2, k3, k4 = paras
    k_IP = np.array([[k0, k1, k2, k3, k4, k4]], dtype=float)
    Pgas = Igas / k_IP#Igas should exclude H2O
    return Pgas
    
def p(I, k, p0):
    #units in 1e-9 mbar and 1e-9 A
    return k*I + p0

def f(x, k):
    return k*x


def readfile(file_name):
    #readfile extracts the data to np arrays. if specific data time needed, use a array.
    Data_labels = {'18': 'H2O', '28': 'N2', '32':'O2', '40': 'Ar', '44': 'CO2', \
         '45': '45CO2', '47':'47CO2', '49':'49CO2', 'TP':'P/mbar', 'H2':'H2'}
    with open(file_name,'r') as f:
        contents = f.read().splitlines()#convert file to list.
        #print(contents[:10])#uncomment to check data
    labels = contents[6].split()#line 7 of the file is where the labels are stored.
    labels = [Data_labels[i] for i in labels]#convert masses and TP to gases and pressure.
    #Source = contents[0]#this line stores the souce filename.
    for i, aline in enumerate(contents):
        contents[i] = aline.split('\t')
    #for contents, from line 9 stores the data.
    aline_len = len(labels)*3
    if aline_len == len(contents[8]):
        print('labels and aline_len match')
    i = -1#usually the last line can be shorter than the rest of the data.
    while len(contents[i]) != aline_len:
        print('there is a bad line:\n', contents[i])
        contents.pop()#remove the bad line
    a = np.array(contents[8:])#list to array
    t_all = np.column_stack([np.asarray(a[:,i], dtype = float) for i in range(1, aline_len, 3)])#stack all t columns
    I_all = np.column_stack([np.asarray(a[:,i], dtype = float) for i in range(2, aline_len, 3)])
    return(t_all, I_all, labels)

def interp(t0, I0, t_new):
    f = interpolate.interp1d(t0, I0, fill_value='extrapolate')
    I_new = f(t_new)
    return I_new

def fit_pI(I, t, p_MS):
    #with I, t for gases, p_MS for pressure measured.
    #the following gives the upper and lower limits for k_PI    
    lower_lims = np.array([0.013, 0.019, 0.0095, 0.00042, 0.0001])
    upper_lims = np.array([0.4, 0.1, 0.1, 0.1, 0.4])
    popt, pcov = curve_fit(I_to_P, I*1e14, p_MS*1e14,bounds = (lower_lims, upper_lims))
    #the reason of using the factor of 1e14 is to make the numericals into "normal range"
    #if the numbers are too small or too large, it won't result a good fit
    #for instance, the square of 1e-8 is already quite small when calculating the fit and data differences.
    plt.plot(t[:,0],p_MS, 'b-', label = 'meaured P/mbar')
    plt.plot(t[:,0],I_to_P(I*1e14, *popt)*1e-14, 'g--', label='fit')
    plt.legend()
    plt.show()
    plt.close()    
    return popt
def log_mag(I, onemin):
    Ilog = np.log2(I)
    dIlog = Ilog[10:]-Ilog[:-10]
    I_relative = I[10:]/np.mean(I[:onemin])
    dIlog_mag = dIlog * I_relative # this is to bring dIlog to close levels
    return dIlog_mag
#dI45log_mag = log_mag(I45, one_min)

def plot_I_changes(I, t, onemin=None, gas = 'unknown', n = 20):
    """

    Parameters
    ----------
    I : np.array (1d)
        current on mass spec reading of a gas species.
    t : np.array
        time, same dimension as I.
    onemin : int, optional
        datapoints within 1 min, can be calculated from t. The default is None.
    gas : str, optional
        the gas that is being analyzed. The default is 'unknown'.
    n : int, optional
        the time range interested to check threshold distinguishing peaks.\
            The default is 20.

    Returns
    -------
    plots for inspection

    """
    if onemin == None:
        onemin = int(round(60/(t[1]-t[0])))
    #onemin is how many data points are within one minute
    first_n_min = n * onemin
    Ilog = np.log2(I)
    dIlog = Ilog[10:]-Ilog[:-10]
    plt.plot(t[:-10], dIlog, 'b-', label = 'dlog2(I) of' + gas)
    show_plot()
    plt.plot(t[:first_n_min], dIlog[:first_n_min], 'b-', \
             label = ' '.join(['dlog2(I) of',gas, str(n),'min']))
    show_plot()
    
    I_relative = I[10:]/np.mean(I[:onemin])
    dIlog_mag = dIlog * I_relative # this is to bring dIlog to close levels
    plt.plot(t[:-10], dIlog_mag, 'g-', label = 'magnified dlog2(I) of' + gas)
    show_plot()

    plt.plot(t[:first_n_min], dIlog_mag[:first_n_min],'g-', \
             label = ' '.join(['magnified dlog2(I) of',gas, str(n),'min']))
    show_plot()
    
    return dIlog_mag

#From: Juan Nunez-Iglesias, Stéfan van der Walt, Harriet Dashnow. “Elegant SciPy.”
def gaussian_kernel(size, sigma):
    """Make a 1D Gaussian kernel of the specified size and standard deviation.

    The size should be an odd number and at least ~6 times greater than sigma
    to ensure sufficient coverage.
    """
    positions = np.arange(size) - size // 2
    kernel_raw = np.exp(-positions**2 / (2 * sigma**2))
    kernel_normalized = kernel_raw / np.sum(kernel_raw)
    return kernel_normalized
def conv_a_sig(sig, size= 25, sigma=4, diff = np.array([1,1,1,0,-1,-1,-1])):
    """
    This function smoothes the change of a noisy signal, to find sharp changes
    """
    smooth_diff = ndi.convolve(gaussian_kernel(size, sigma), diff)
    sdsig = ndi.convolve(sig, smooth_diff)#smoothed dsig
    plt.plot(sdsig)
    return sdsig

def find_rises(dC45, dImin = [1,1,1,10,10,10,100]):
    """
    an alternative away to find time points where bicarbonate were added
    gives imaxes which is the index numbers indentified.

    Parameters
    ----------
    dC45 : numpy.array
        '45CO2' signal, convolved.
        an example given:
        #dC45 = conv_a_sig(I_df['45CO2'], size=63, diff=np.array([1]*30+[0]+[-1]*30))
        #dC45 = dC45*1e12
    dImin : list, optional
        a list of values for identifying rises, plot dC45 upto different index
        number will help to determine the proper dImin to use. 
        The default is [1,1,1,10,10,10,100].

    Returns
    -------
    imaxes, a list of indexes for identifying points.
    NOTE: one may need to re-slice and adjust the points becasue the list returned
    usually are on the slope, not at the base

    """
    imaxes = []
    idx_start =0
    onemin = 17
    for j in range(len(dImin)):
        idx = np.where(dC45[idx_start:]> dImin[j])[0][0] + idx_start
        idxmax = np.where(dC45==np.max(dC45[idx:idx+onemin*2]))[0][0]
        idx_start = idxmax + onemin*2
        #index increment for slicing np.array dI10   
        """
        the next idx for stable pressure should come after ~3 minutes,
        2 min at least to find the next. 
        """    
        imaxes.append(idxmax)
        dC45[:idxmax+onemin*2]=0
        plt.plot(dC45, label = 'after point'+str(j+1))
        show_plot()
        # plt.plot(I_df['45CO2'].iloc[idxmax-15:idxmax+15])
        # show_plot()
    i_end = np.argwhere(dC45==dC45.min())[0][0]
    imaxes.append(i_end)
    return imaxes


def calibrate_CO2(I45, tp, onemin, thresholds = [1,10,100]):
    dI10 = log_mag(I45, onemin)
    imaxes = []

    #plot_I_changes(I45, t_p, one_min, '45CO2')
    for i in range(3):
        dImin = thresholds[i]
        #the dImin is used to identify peaks/ time points of adding NaHCO3
        idx_start = 0
        if i == 2:#when i==2, do the following once
            idx = np.where(dI10[idx_start:]> dImin)[0][0] + idx_start
            idxmax = np.where(dI10==np.max(dI10[idx:idx+onemin]))[0][0]
            idx_start = idxmax + onemin*2
            #index increment for slicing np.array dI10   
            """
            the next idx for stable pressure should come after ~3 minutes,
            2 min at least to find the next. 
            """    
            imaxes.append(idxmax)
            dI10[:idxmax+17]=0#make the identified spikes as 0, 17 points is
            #half a min in most cases and no longer than one min
        else:
            for j in range(3):
                idx = np.where(dI10[idx_start:]> dImin)[0][0] + idx_start
                idxmax = np.where(dI10==np.max(dI10[idx:idx+onemin]))[0][0]
                idx_start = idxmax + onemin*2
                #index increment for slicing np.array dI10   
                """
                the next idx for stable pressure should come after ~3 minutes,
                2 min at least to find the next. 
                """    
                imaxes.append(idxmax)
                dI10[:idxmax+17]=0

    i_end = np.where(dI10 == np.min(dI10))[0][0] - 2   
    imaxes.append(i_end)
    #imaxes.append(-1)#use this if the end is equilibrium state
    x = tp[imaxes]
    y = np.log10(I45[imaxes])
    plt.plot(tp, np.log10(I45), 'b-', label = 'data of log10(I_45CO2)')
    plt.plot(x, y, 'ro', label = 'identified points')
    show_plot()
    yI = []
    for i in imaxes:
        yI.append(np.mean(I45[i-13:i-3]))
    return (imaxes, yI)

def fit_I_C(xCO2, c):
    
    popt3, pcov3 = curve_fit(f, xCO2*1e12, c*1e3)#pA vs uM
    plt.plot(xCO2*1e12, c*1e3, 'bo', label = 'dc vs dI data')
    plt.plot(xCO2*1e12, f(xCO2*1e12, *popt3), 'g--', label = 'c = %5.2f * I' %tuple(popt3))
    plt.xlabel('I /pA')
    plt.ylabel('$^{13}CO_2$ /µM')
    show_plot()
    return popt3[0]

def calc_f45(I44, I45, tp, i_maxes):
    #calculate the fraction of 45CO2 in the 13C labeled NaHCO3
    plt.plot(tp, I44, label = '44CO2')
    show_plot()
    
    I0 = np.mean(I44[:30])    
    dI10 = 10*(I44[10:]-I44[:-10])/I0
    plt.plot(tp[:-10], dI10, 'bo', markersize = 2, mfc = 'None', label = 'dI10/I0 of 44CO2')
    show_plot()
    #i = np.where(dI10>1)[0][0]
    #i = np.where(np.max(dI10) == dI10)[0][0]
    """the1 value should be changed based on the observation of dI10 curves"""
    i_end = i_maxes[-1]
    i_0 = i_maxes[-2]
    dI44 = np.mean(I44[i_end-20:i_end])-np.mean(I44[i_0-10:i_0])
    dI45 = np.mean(I45[i_end-20:i_end])-np.mean(I45[i_0-10:i_0])
    f45 = dI45/(dI44+dI45)#0.9856
    return f45

def file_to_DF(filename, plot = True):
    """
    read a file and return a DataFrame of current on mass spec
    interpolated to synchronize the datapoints
    Parameters
    ----------
    filename : str
        MIMS datafile.
    plot : Bool, optional
        whether or not plot the CO2 O2 data. The default is True.

    Returns
    -------
    I_df : Pandas.DataFrame
        DataFrame using mass 45 as time reference.

    """
    t_all, I_all, labels = readfile(filename)#assign t, I, labels
    
    index_45 = labels.index('45CO2')
    t_45 = t_all[:, index_45]#use the 45CO2 as the time for interpolation and calculations    
    #one can change which time index to use, here using mass 45

    I_gas = {}    
    for i, a_gas in enumerate(labels[:-1]):
        #store current in to I_gas dictionary, interpolated using linear method
        I_gas[a_gas] = interp(t_all[:,i], I_all[:,i], t_45)
    #store interpolated data in to Pandas.DataFrame
    I_df = pd.DataFrame(I_gas, index = t_45)
    
    ####plot the data to see the change of CO2 signal by O2####
    if plot:
        if '47CO2' in I_df:
            CO2s = ['45CO2', '47CO2','49CO2']
        else:
            CO2s =['45CO2']
        I_df.plot(y=CO2s)
        ax = I_df['O2'].plot(secondary_y = True, color ='b', marker = 'o', markersize = 1)
        ax.set_ylabel('$O_2$', color = 'b')
        plt.show()
        
        I_df.plot(y=['CO2'])
        ax = I_df['O2'].plot(secondary_y = True, color ='b', marker = 'o', markersize = 1)
        ax.set_ylabel('$O_2$', color = 'b')
        plt.show()
        plt.close()
    return I_df

def O2_sol(t, S):
    """
    Calculate O2 solubility in water/seawater
    according to 
    eq8 in Garcia and Gordon 1992 Limol. Oceanogr. 37(6) 1307-1312
    Oxygen solubility in seawater: Better fitting equations
    NOTE: tF≤t≤40°C, 0≤S≤42 [typo in original article]
    &
    eq2 in Mortimer 1981
    The oxygen content of air-saturated fresh waters
    over ranges of temperature and atmospheric
    pressure of limnological interest
    Parameters
    ----------
    t : TYPE
        Temperature in °C.
    S : TYPE
        Salinity in per mil.

    Returns
    -------
    Air saturated Oxygen concentration in µmol/kg, if S=0, µM.

    """
    A0 = 5.80818
    A1 = 3.20684
    A2 = 4.11890
    A3 = 4.93845
    A4 = 1.01567
    A5 = 1.41575
    
    B0 = -7.01211e-3
    B1 = -7.25958e-3
    B2 = -7.93334e-3
    B3 = -5.54491e-3
    
    C0 = -1.32412e-7
    #the following calculate O2 concentration in water/seawater
    if S!=0:
        Ts = np.log((298.15-t)/(273.15+t))
        lnC = A0 + A1*Ts + A2*Ts**2 + A3*Ts**2 + A3*Ts**3\
            + A4*Ts**4 + A5*Ts**5\
                +S*(B0 + B1*Ts + B2*Ts**2 + B3*Ts**3) + C0*S**2
        O2_solubility = np.exp(lnC)
    else:
        T = t+273.15
        lnC= -139.34410 + (1.575701e5 /T)- (6.642308e7 /T**2)\
            + (1.243800e10/T**3 )- (8.621949e11 /T**4 )
        O2_solubility = np.exp(lnC)
        O2_solubility = O2_solubility*1000/32#convert mg/L to µM
    
    return O2_solubility


def I_to_c(I_df, Imins):
    c = pd.DataFrame()
    kO2 = Imins['kO2']
    kIC = Imins['kIC']
    f13C= Imins['f13C']
    if 'O2min' in Imins:
        c['O2'] =(I_df['O2'] -Imins['O2min'])*kO2
    else:
        c['O2'] =(I_df['O2'] -Imins['O2'])*kO2
    I_df['netCO2'] = I_df['CO2']-\
        (I_df['O2'] -Imins['O2'])*Imins['d44dO2']- Imins['CO2']
    I_df['net45CO2'] = I_df['45CO2']-\
        (I_df['O2'] -Imins['O2'])*Imins['d45dO2']- Imins['45CO2']
    if '47CO2' in I_df:
        I_df['net47CO2'] = I_df['47CO2'] - Imins['47CO2']
        I_df['net49CO2'] = I_df['49CO2'] - Imins['49CO2']
    #the following assuming background CO2 at 0    
    # popt3, pcov3 = curve_fit(f, I_df[I_df.index<2000]['O2'], \
    #                          I_df[I_df.index<2000]['netCO2'])
    # # index 2000s is where labeled bicarbonate added.
    # I_df['netCO2'] = I_df['netCO2'] - (I_df['O2'] -Imins['O2'])* popt3
    # I_df['net45CO2'] = I_df['45CO2']-\
    #     (I_df['O2'] -Imins['O2'])*Imins['d45dO2']
    # #due to low impact for 45CO2, one can ignore the impact?
    if '47CO2' in I_df:
        CO2s = ['CO2', '45CO2', '47CO2', '49CO2']
    else:
        CO2s = ['CO2', '45CO2']
    for aCO2 in CO2s:
        c[aCO2] = I_df['net'+aCO2] *kIC*f13C*1e12
    
    # I_df.plot(y=['netCO2', 'net45CO2', '47CO2', 'net49CO2'])
    # ax = I_df['O2'].plot(secondary_y = True, color ='b', marker = 'o', markersize = 1)
    # ax.set_ylabel('$O_2$', color = 'b')
    # plt.show()
    
    c.plot(y=CO2s)
    plt.ylabel('$CO_2$/µM')
    ax = c['O2'].plot(secondary_y = True, color ='b', marker = 'o', markersize = 1)
    ax.set_ylabel('$O_2$', color = 'b')
    plt.show()    
    return c

def plot_df(df, CO2 = 'CO2',O2=True, total= False):
    """
    Plot a DataFrame where it has 13C labeled CO2
    ----------
    df : DataFrame
       DataFrame of CO2, O2 concentrations.
    O2: bool; whether or not plot O2
    -------
    a plot of CO2 (and O2 if O2 == True)

    """
    df.plot(y=CO2, style='o', markersize = 2)
    plt.ylabel('$CO_2$/µM')
    if O2:
        ax = df['O2'].plot(secondary_y = True, color ='b', marker = 'o', markersize = 1)
        ax.set_ylabel('$O_2$', color = 'b')
    if total:
        total_co2=df[['45CO2','47CO2','49CO2']].sum(axis=1)
        ax = total_co2.plot(color = 'k', marker = 'o', markersize = 2)
        ax.legend(['45CO2','47CO2','49CO2','total_13CO2'], loc= 'best')
    plt.rcParams["figure.dpi"] = 300
    plt.legend()
    plt.show()
    
def fit_O2_Ar(df, length):
    """
    fit O2 and Ar signal to a linear line for O2 correction

    Parameters
    ----------
    df : DataFrame
        DESCRIPTION.
    length : int
        The index of DataFrame of how far the dataset is fitting.

    Returns
    -------
    slope, intercept and r_square.

    """
    m32 = df[df.index<length]['O2'].to_numpy()
    m40 = df[df.index<length]['Ar'].to_numpy()
    slope, intercept, r_value, p_value, std_err=stats.linregress(m40,m32)
    r_square = r_value**2
    plt.plot(m40, m32, 'bo', markersize = 2, label ='data')
    plt.plot(m40, m40*slope+intercept, 'g--', label = 'fit: slope = %5.5f' %slope)
    plt.xlabel('Ar_m40')
    plt.ylabel('$O_2$_m32')
    show_plot()
    return (slope, intercept, r_square)

#### calculate minimum CO2 currents and dCO2/dO2 ====####
filename = 'CO2_wo_47_49_CO2_calib.asc'
#This can also be the file that stores data from adding Na2S2O4 to 0.1 N NaOH
df1 = file_to_DF(filename, plot = True)#df1 for NaOH, Na2S2O4
df1[df1.index>2500].plot(y='CO2')

df_temp = df1[df1.index>2500]
plot_df(df_temp)

Imins = df1[df1.index>3200][['O2','CO2','45CO2']].mean()
Imins = df1[df1.index>3200][['O2','CO2','45CO2']].mean()
#there is no mass 47 and 49 in this file, but mass 47 and 49 \
#are expected to be used later. mass 45, 47, 49 are treated as zero \
#before adding 18O-labeled NaHCO3
Imins['47CO2']=Imins['45CO2']#this is for later calculation
Imins['49CO2']=Imins['45CO2']#this is for later calculation
df_temp.iloc[:100].plot()
I0s = df_temp.iloc[:100].mean()

dIs = I0s-Imins
Imins['d45dO2']=dIs['45CO2']/dIs['O2']
Imins['d44dO2']=dIs['CO2']/dIs['O2']
#Imins.to_csv(filename[:-4]+' mins.csv')
#unquote the above to save Imins

#### calibrate O2,  process H2O followed by Na2S2O4 ====####
f2 = 'CO2_wo_47_49_O2_calib.asc'
df2 = file_to_DF(f2)#df2 for H2O Na2S2O4
O2min = df2['O2'][df2.index>800].mean()
df2.iloc[100:200].plot()
O2max = df2.iloc[100:200]['O2'].mean()
dO2 = O2max-O2min
O2_sat = O2_sol(3,0)# this will be different at different temperature
kO2 = O2_sat/dO2

Imins['O2max'] = O2max
Imins['kO2'] = kO2
Imins['O2min'] = O2min
#Imins.to_csv(filename[:-14]+' mins.csv')
#unquote the above to save Imins

#### calibrate CO2 concentration vs signal ====#### 
f_calib = 'CO2_wo_47_49_CO2_calib.asc'
I_df = file_to_DF(f_calib)

t45 = I_df.index
I45 = I_df['45CO2'].to_numpy()
one_min = int(round(60/(t45[1]-t45[0])))
# plot_I_changes(I45, t45, one_min, '45CO2', 15)
# plot_I_changes(I45, t45, one_min, '45CO2', 35)
# plot_I_changes(I45, t_p, one_min, '45CO2', 30)
#The above plots help determine what thresholds should be used to identify


dC45 = conv_a_sig(I_df['45CO2'], size=63, diff=np.array([1]*30+[0]+[-1]*30))
# this helps the identification of time points adding NaH13CO3
dC45 = dC45*1e12
plt.plot(dC45[:1000])#show first 1000 data points of convolved signal.
plt.plot(dC45[:500])#helps to see if more than 7 steps
show_plot()

thresholds =[1.0, 1.5, 1.0, 15, 12, 12, 150]#adjust the values and rerun from dC45
imaxes = find_rises(dC45, thresholds)

# imaxes = imaxes[:4]+imaxes[5:]
imaxes = np.array(imaxes)-18
# 8-point after the base is common, one can change the number to see if the 
# indentification is accurate by plot the following data
# imaxes[3] = imaxes[3]-90#this is to change individual data points
# imaxes[1] = imaxes[1]-60
x = I_df.index[imaxes]
y = np.log10(I45[imaxes])
plt.plot(t45, np.log10(I45), 'b-', label = 'data of log10(I_45CO2)')
plt.plot(x, y, 'ro', label = 'identified points')
show_plot()
ICO2 = []
for i in imaxes:
    ICO2.append(np.mean(I45[i-10:i]))

"""" alternatively use the following lines to get imaxes and ICO2"""

# thresholds = [5,10,200]#see the plot of dI45log_mag for thresholds
# imaxes, ICO2 = calibrate_CO2(I45, t45, one_min, thresholds)

c = np.array([0, 3/1203, 6/1206, 9/1209, 39/1212, 69/1215, 99/1218,\
              399/1221])#unit mM
xCO2 = np.array(ICO2)-ICO2[0]#xCO2 was a list of numpy.float, but now it can be\
# converted to array because xCO2[0] is numpy.float? current in picoA
kIC = fit_I_C(xCO2, c)

# xCO2 = xCO2[:7]
# c = c[:7]
kIC = fit_I_C(xCO2[:-3], c[:-3])#calibrate the I-c relation of CO2, with plot
# c[CO2] = kIC*I[CO2]*1e12, units in µM for c, A for I

I44 = I_df['CO2'].to_numpy()
f13C = calc_f45(I44, I45, t45, imaxes)
Imins['kIC'] = kIC 
Imins['f13C'] = f13C  
Imins.to_csv(filename[:-14]+' mins_paras.csv')

#### convert I to c for exp file ####
f_exp = 'Fcyl_3C_S40_exp01_wAZ.asc'
df_exp = file_to_DF(f_exp)

O2_Ar=(df_exp['O2']/df_exp['Ar']).iloc[:500].mean()#initial O2/Ar signal
(df_exp['O2']-df_exp['Ar']*O2_Ar).plot()

slp, intc, r2 = fit_O2_Ar(df_exp, 500)

#correct O2 based on Ar signal
df_exp['O2_Ar']=df_exp['O2']-df_exp['Ar']*O2_Ar + df_exp['O2'].iloc[:-30].mean()

#df_exp['O2_Ar']=df_exp['O2']-df_exp['Ar']*slp-intc + df_exp['O2'].iloc[-30:].mean()


cs = I_to_c(df_exp, Imins)#concentrations in DataFrame
cs['O2_Ar'] = (df_exp['O2_Ar'] -Imins['O2'])*kO2
cs.to_csv(f_exp[:-14]+'_exp01_c.csv')

#f13C= 0.98158312    
#### convert exp02 I to c####
f_exp2 = 'Fcyl_3C_S40_exp02_noAZ.asc'
df_exp2 = file_to_DF(f_exp2)

O2_Ar=(df_exp2['O2']/df_exp2['Ar']).iloc[:500].mean()#initial O2/Ar signal
(df_exp2['O2']-df_exp2['Ar']*O2_Ar).plot()

slp, intc, r2 = fit_O2_Ar(df_exp2, 500)

#correct O2 based on Ar signal

df_exp2['O2_Ar']=df_exp2['O2']-df_exp2['Ar']*O2_Ar + df_exp2['O2'].iloc[:-30].mean()
#df_exp2['O2_Ar']=df_exp2['O2']-df_exp2['Ar']*slp - intc + df_exp2['O2'].iloc[-30:].mean()


cs2 = I_to_c(df_exp2, Imins)#concentrations in DataFrame
cs2['O2_Ar'] = (df_exp2['O2_Ar'] -Imins['O2'])*kO2
cs2.to_csv(f_exp2[:-14]+'_exp02_c.csv')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
