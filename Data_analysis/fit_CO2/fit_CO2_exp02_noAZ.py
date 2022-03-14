#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 10:33:29 2020
Fit data to CO2 related kinetics
@author: LIMeng, limco2(AT)uw.edu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages#save multi fig to pdf
from scipy.optimize import curve_fit
from scipy import stats
from scipy.integrate import solve_ivp
from scipy.special import comb
from scipy import ndimage as ndi
#from scipy.interpolate import splev, splrep
from scipy.signal import savgol_filter

def show_plot():
    plt.legend()
    plt.show()
    plt.close()
def plot_df(df, O2=True, total= False):
    """
    Plot a DataFrame where it has 13C labeled CO2
    ----------
    df : DataFrame
       DataFrame of CO2, O2 concentrations.
    O2: bool; whether or not plot O2
    -------
    a plot of CO2 (and O2 if O2 == True)

    """
    if 'O2_Ar' in df:
        oxygen = 'O2_Ar'
    else:
        oxygen = 'O2'
    df.plot(y=['45CO2', '47CO2', '49CO2'], style='o', markersize = 2)
    plt.ylabel('$CO_2$/µM')
    if O2:
        ax = df[oxygen].plot(secondary_y = True, color ='b', marker = 'o', markersize = 1)
        ax.set_ylabel('$O_2$', color = 'b')
    if total:
        total_co2=df[['CO2','45CO2','47CO2','49CO2']].sum(axis=1)
        ax = total_co2.plot(color = 'k', marker = 'o', markersize = 2)
        ax.legend(['45CO2','47CO2','49CO2','total_CO2'], loc= 'best')
    plt.rcParams["figure.dpi"] = 300
    plt.legend()
    plt.show()

def plot_totalC(df, O2 = True):
    if 'O2_Ar' in df:
        oxygen = 'O2_Ar'
    else:
        oxygen = 'O2'
    df[['CO2','45CO2','47CO2','49CO2']].sum(axis = 1).plot(style='-', markersize = 2)
    plt.ylabel('$total [CO_2]$/µM')
    if O2:
        ax = df[oxygen].plot(secondary_y = True, color ='b', marker = 'o', markersize = 1)
        ax.set_ylabel('$O_2$', color = 'b')
    plt.legend()
    plt.show()
    plt.close()    

def pK(TC, S):
    T = TC + 273.15
    pK1= -43.6977- 0.0129037*S + 1.364e-4*S**2 +2885.378/T+ 7.045159*np.log(T)
    pK2= -452.0940+ 13.142162*S -8.101e-4*S**2 +21263.61/T+ 68.483143*np.log(T)\
        +(-581.4428*S +0.259601*S**2)/T - 1.967035*S * np.log(T)
    return (pK1, pK2)

def visc(Tc, S, V = False):
    """
    Calculate the viscosity of water or seawater
    Reference:
        Sharqawy et al. 2009
        Thermophysical properties of seawater: a review of \
            existing correlations and data
        https://doi.org/10.5004/dwt.2010.1079
        Eq 22 and Eq 23
        viscosity unit in kg/m/s
    Parameters
    ----------
    Tc : scalar
        Temperature in °C.
    S : scalar
        Salinity in per mil (g/kg).
    V:  Bool
        whether or not return calculated viscosity, default is False

    Returns
    -------
    Viscosity or Seawater to water ratio.

    """
    Sref = S/1000 #reference salinity is written in kg/kg
    A = 1.541 + 1.998e-2 * Tc - 9.52e-5 * Tc**2
    B = 7.974 - 7.561e-2 * Tc + 4.724e-4 *Tc**2
    Sw_w = 1 + A*Sref + B*Sref**2#seawater to water ratio

    if V:        
        V_water = 4.2844e-5 + 1/( 0.157* (Tc+ 64.993)**2 -91.296)
        #water viscosity
        V_sw = V_water * Sw_w # seawater viscosity
        Sw_w = (Sw_w, V_sw)
        print("given viscosity is in kg/m/s, which is 1000 times of centipoise")
    
    return Sw_w

def Di(Tc, S=0, Ci=None):
    """
    Calculate Diffussion coefficient of CO2, HCO3-, CO3--, according to 
    R.E. Zeebe / Geochimica et Cosmochimica Acta 75 (2011) 2483–2498
    Original equation covers Tc in [0,100]

    Parameters
    ----------
    Tc :scalar
        temperature in  °C.
    S: scalar
        salinity in per mil (or g/kg)
    Ci : string, optional, default is None
        CO2, HCO3-, CO3--

    Returns
    -------
    A value of diffussion coefficient of a caron species 
    or a dictionary if Ci is not specified.

    """
    Ci_list = ['CO2', 'HCO3-','CO3--']
    D0i = np.array([14.6836, 7.0158, 5.4468])*1e-5
    #unit in cm**2/s, or change 1e-5 to 1e-9 for m**2/s
    Ti  = np.array([217.2056, 204.0282, 210.2646])
    gamma_i = np.array([1.997, 2.3942, 2.1929])
    T = Tc + 273.15
    D = D0i * ((T/Ti)-1) **gamma_i
    if Ci in Ci_list:
        i = Ci_list.index(Ci)
        D = D[i]
    if S!=0:
        D=D/visc(Tc, S)

    return D


def f_passive(t, y, kuf, kur):
    #the simple differential equations for HCO3- --> CO2 kinetics
    """    
    CO2   --kuf--> HCO3-
    HCO3- --kur--> CO2
    #note that above kinetics upholds overtime only when pH is constant
    
    Parameters
    ----------
    y : tuple
        store all changing parameters: CO2, HCO3- in 
        different isotopic compositions
        Bic, CO3: bicarbonate, carbonate conctrations
    t : float
        time.
    kuf, kur : float
        forward and reverse rxn kinetic constants.
        u for uncatalyzed
    Returns
    -------
    dydt list of variables

    """
    global G, H
    C45, C47, C49, B62, B64, B66, B68 = y
    CO2 = np.array([C45, C47, C49])
    #CO2 is a (3, 1) array, mass 45, 47, 49; or matrix if in matlab
    Bic = np.array([B62, B64, B66, B68])
    #Bic is a (4, 1) array, mass [[m62], [m64], [m66], [m68]]
    dCO2 = -kuf * CO2 + kur * (H @ Bic)
    dBic =  kuf * (G @ CO2) - kur * Bic
    dC45, dC47, dC49 = dCO2
    dB62, dB64, dB66, dB68 = dBic
    dydt = [dC45, dC47, dC49, dB62, dB64, dB66, dB68]
    return dydt

def sim_passive(T_span, kuf, kur, DIC_t = 2072.6, f13C = 0.985, f18O = 0.97, T_eval=None):
    """
    simulate passive dehydration of 18O-HCO3-
    Parameters
    ----------
    T_span : a tuple of two time float
        start and end time for the simulation.
    kuf : float
        forward rxn k, CO2 --kuf-->HCO3-.
    kur : float
        reverse rxn k, HCO3- --kur-->CO2.
    DIC_t : float/int, optional
        DIC added as total initial concentration, in µM. 
        The default is 2500/1.2062 (µM).
    f13C: the fraction of 13C in total C in 13C-NaHCO3.
        The default is 0.985, of our chemical stock.
    T_eval: 1d array of time in seconds, used for in solve_ivp
        The default is None

    Returns
    -------
    Simulated isotopic CO2 concentrations.

    """
    y0 = np.zeros((7,))
    #0, 1, 2 are mass 45, 47, 49 respectively
    # 3-->6, are HCO3- with mass 62, 64, 66, 68
    for i in range(4):
        y0[3+i] = DIC_t * comb(3, i) * (f18O**i) * ((1-f18O)**(3-i)) * f13C
        """what is f18O? 0.8835 as max tau value 18O% in CO2?"""
        #comb, is scipy.special.comb(), calculate number of combinations
    
    sol = solve_ivp(f_passive,T_span,y0, args = (kuf, kur), method = "BDF", \
                t_eval=T_eval, rtol=1e-6, atol=1e-10, max_step=5)
    #plt.plot(sol.t,sol.y[0,:]) #uncomment to see the each tbrk.
    return sol

def fit_passive(passive_df, F18O = 0.8835, F13C = 0.985):
    def func(t_array, kuf0, kur0):
        s = sim_passive([passive_df.index[0], passive_df.index[-1]],\
                        kuf0, kur0, f13C = F13C, f18O = F18O, T_eval = t_array)
        co2s = s.y[:3, :]
        return co2s.flatten()
    low_lims = np.array([1.0e-4, 1.0e-6])
    up_lims = np.array([1.0, 1.0e-2])
    init_guess = np.array([4.0e-2, 4.0e-4])
    ydata = passive_df[['45CO2','47CO2','49CO2']].T.to_numpy().flatten()
    #sim_passive return (3, N) array
    xdata = passive_df.index.to_numpy()#time array
    popt1, pcov1 = curve_fit(func, xdata, ydata, \
                              p0 =init_guess,bounds = (low_lims, up_lims))
        
    passive_df.plot(y=['45CO2','47CO2','49CO2'], style ='o', markersize = 2)
    co2s = func(xdata,*popt1).reshape(3,-1)
    labels = ['C45', 'C47', 'C49']
    for i in range(3):
        plt.plot(xdata,co2s[i,:], '--', label = labels[i]+'fit')
    show_plot()
    return popt1

def fit_passive_3(passive_df, F13C = 0.985):
    def func(t_array, kuf0, K_u, DIC, F18O):
        kur0 = kuf0/K_u
        s = sim_passive([passive_df.index[0], passive_df.index[-1]],\
                        kuf0, kur0, DIC_t = DIC, f13C = F13C, f18O = F18O, T_eval = t_array)
        co2s = s.y[:3, :]
        return co2s.flatten()
    low_lims = np.array([4.20e-3, 104.7, 1900.0, 0.75])
    #limits of kuf0, kur0, DIC in µM, F18O
    up_lims = np.array([5.26e-3, 127.9, 2100.0, 0.984])
    init_guess = np.array([4.78e-3, 116.3, 2073.0, 0.983])
    ydata = passive_df[['45CO2','47CO2','49CO2']].T.to_numpy().flatten()
    #sim_passive return (3, N) array
    xdata = passive_df.index.to_numpy()#time array
    popt1, pcov1 = curve_fit(func, xdata, ydata, \
                              p0 =init_guess,bounds = (low_lims, up_lims))
        
    passive_df.plot(y=['45CO2','47CO2','49CO2'], style ='o', markersize = 2)
    co2s = func(xdata,*popt1).reshape(3,-1)
    labels = ['C45', 'C47', 'C49']
    for i in range(3):
        plt.plot(xdata,co2s[i,:], '--', label = labels[i]+'fit')
    show_plot()
    return popt1

def f_dark_SF(t, y, ksf, kcf, KCYTO):
    """
    Boundary Layer Model, dark
    Eqs.1-6 from Hopkinson et al. June 2013, Vol. 162, pp. 1142–1152
    Note that units are in µM, mL. fc, fb in mL/s per cell
    All k values are in s-1 unit, 
    the ksf in unit of cm3/s (as in reference) needs another step of calculation
    """
    #C for CO2, B for bicarbonate, \
    #e,s,i for extracellular, surface, and intracellular
    C45e, C47e, C49e, B62e, B64e, B66e, B68e,\
    C45s, C47s, C49s, B62s, B64s, B66s, B68s,\
    C45i, C47i, C49i, B62i, B64i, B66i, B68i = y
    global G, H, k1, k2, N, Ve, Vs, Vc, fc_BL, fb_BL, fc_M, fb_M
    ksr = k2*ksf/k1#KSF is the surface equilibrium constant
    kcr = ksf/KCYTO#KCYTO is the cytosolic equilibrium constant
    Ce = np.array([C45e, C47e, C49e])
    Be = np.array([B62e, B64e, B66e, B68e])
    Cs = np.array([C45s, C47s, C49s])
    Bs = np.array([B62s, B64s, B66s, B68s])
    Ci = np.array([C45i, C47i, C49i])
    Bi = np.array([B62i, B64i, B66i, B68i])
    #Eqs.1-6
    dCe = -k1 * Ce + k2* (H @ Be) + fc_BL * N *(Cs - Ce)/Ve
    dBe = k1 * (G @ Ce) - k2 * Be + fb_BL * N *(Bs - Be)/Ve
    # dCs = -ksf* Cs/(Cs+Kcs) + ksr*(H @ (Bs/(Bs+Kbs))) + (1/Vs)* (fc_BL* (Ce - Cs) + fc_M*(Ci - Cs))
    # dBs = ksf* (G @ (Cs/(Cs+Kcs))) - ksr* Bs/(Bs+Kbs) + (1/Vs)* (fb_BL* (Be - Bs) + fb_M*(Bi - Bs))
    dCs = -ksf* Cs + ksr*(H @ Bs) + (1/Vs)* (fc_BL* (Ce - Cs) + fc_M*(Ci - Cs))
    dBs = ksf* (G @ Cs) - ksr* Bs + (1/Vs)* (fb_BL* (Be - Bs) + fb_M*(Bi - Bs))
    #here the ksf and ksr do not include Vs term, so to get the ksf value in reference, ksf_ref = ksf * Vs
    #The Vs in Hopkinson 2013 Plant Phys should be Vs/per_cell * N, as total surface volume.
    
    dCi = -kcf* Ci + kcr *(H @ Bi)+ fc_M * 1 *(Cs - Ci)/Vc#[S3]
    dBi = kcf* (G @ Ci) - kcr* Bi + fb_M * 1 *(Bs - Bi)/Vc#[S4]
    #NOTE that the equation Vi in the original paper must refer to total cytosolic
    #volume of cells, i.e. Vi = N*Vc
    
    # each individual CO2 and HCO3- species:
    dC45e, dC47e, dC49e = dCe
    dB62e, dB64e, dB66e, dB68e = dBe
    dC45s, dC47s, dC49s = dCs
    dB62s, dB64s, dB66s, dB68s = dBs
    dC45i, dC47i, dC49i = dCi
    dB62i, dB64i, dB66i, dB68i = dBi
    dydt = np.concatenate((dCe, dBe, dCs, dBs, dCi, dBi)).tolist()
    #dydt = [dC45e, dC47e, dC49e, dB62e, dB64e, dB66e, dB68e,\
    #        dC45s, dC47s, dC49s, dB62s, dB64s, dB66s, dB68s,\
    #        dC45i, dC47i, dC49i, dB62i, dB64i, dB66i, dB68i]
    return dydt

def sim_dark_SF(T_span, Y0, ksf, kcf, KCYTO, T_eval = None):
    sol_dark = solve_ivp(f_dark_SF,T_span,Y0, method='BDF',\
                         args = (ksf, kcf, KCYTO), t_eval=T_eval)
    return sol_dark

def fit_dark_SF(dark_df, Y0):
    """ the Y0 has 7 x 3 = 21 parameters, see f_dark_SF"""
    def f_SF(t_array, ksf, kcf, KCyto):
        s = sim_dark_SF([dark_df.index[0], dark_df.index[-1]], Y0, \
                        ksf, kcf, KCyto,  T_eval = t_array)
        co2s = s.y[:3, :]
        #print(s.y.shape)
        return co2s.flatten()
    #to get ksf, ksr
    low_lims = np.array([25, 20, 9.21])
    up_lims = np.array([500, 170, 92.1 ])
    init_guess = np.array([164, 85.7, 68.6])
    ydata = dark_df[['45CO2','47CO2','49CO2']].T.to_numpy().flatten()
    #print(ydata.shape, type(ydata))
    #to fit ydata, (3,N) array cannot be calculated, so fatten
    xdata = dark_df.index.to_numpy()#time array
    #print(xdata.shape, type(xdata))
    #j = f_dark(xdata, *init_guess)
    #print(j.shape, type(j))
    popt2, pcov2 = curve_fit(f_SF, xdata, ydata, \
                              p0 =init_guess,bounds = (low_lims, up_lims),\
                                  loss = 'soft_l1')
        
    dark_df.plot(y=['45CO2','47CO2','49CO2'], style ='o', markersize = 2)
    CO2s = f_SF(xdata,*popt2).reshape(3,-1)
    labels = ['C45', 'C47', 'C49']
    for i in range(3):
        plt.plot(xdata,CO2s[i,:], '--', label = labels[i]+'fit')
    show_plot()
    return popt2    

def f_dark_1C(t, y, kcf, kcr, fc, fb, eCA):
    """
    One-Compartment Model, dark
    Eqs.S1-S4 from Hopkinson et al. 10.1073/pnas.1018062108
    Note that units are in µM, mL. fc, fb in mL/s per cell
    
    eCA is the factor that accelerate k1, k2
    """
    #C for CO2, B for bicarbonate, e, i for extracellular and intracellular
    C45e, C47e, C49e, B62e, B64e, B66e, B68e,\
    C45i, C47i, C49i, B62i, B64i, B66i, B68i = y
    global G, H, k1, k2, N, Ve, Vc
    Ce = np.array([C45e, C47e, C49e])
    Be = np.array([B62e, B64e, B66e, B68e])
    Ci = np.array([C45i, C47i, C49i])
    Bi = np.array([B62i, B64i, B66i, B68i])
    #Eqs.S1-S4
    dCe = -k1 * Ce * eCA + k2* (H @ Be) * eCA + fc * N *(Ci - Ce)/Ve
    dBe = k1 * (G @ Ce) * eCA - k2 * Be * eCA + fb * N *(Bi - Be)/Ve
    dCi = -kcf* Ci + kcr *(H @ Bi)+ fc * 1 *(Ce - Ci)/Vc#[S3]
    dBi = kcf* (G @ Ci) - kcr* Bi + fb * 1 *(Be - Bi)/Vc#[S4]
    #NOTE that the equation Vi in the original paper must refer to total cytosolic
    #volume of cells, i.e. Vi = N*Vc
    
    # each individual CO2 and HCO3- species:
    dC45e, dC47e, dC49e = dCe
    dB62e, dB64e, dB66e, dB68e = dBe
    dC45i, dC47i, dC49i = dCi
    dB62i, dB64i, dB66i, dB68i = dBi
    dydt = np.concatenate((dCe, dBe, dCi, dBi)).tolist()
    #dydt = [dC45e, dC47e, dC49e, dB62e, dB64e, dB66e, dB68e,\
    #        dC45i, dC47i, dC49i, dB62i, dB64i, dB66i, dB68i]
    return dydt

def sim_dark_1C(T_span, Y0, kcf, kcr, fc, fb, eCA, T_eval = None):
    sol_dark = solve_ivp(f_dark_1C,T_span,Y0, method='BDF',\
                         args = (kcf, kcr, fc, fb, eCA), t_eval=T_eval)
    return sol_dark
    
def fit_dark_1C(dark_df, Y0):
    """ the Y0 has 7 x 2 = 14 parameters, see f_dark_1C"""
    def f_dark(t_array, kcf, kcr, fc, fb, eCA):
        s = sim_dark_1C([dark_df.index[0], dark_df.index[-1]], Y0, \
                        kcf, kcr, fc, fb, eCA, T_eval = t_array)
        co2s = s.y[:3, :]
        #print(s.y.shape)
        return co2s.flatten()
    #to get kcf, kcr, fc, fb, and eCA
    low_lims = np.array([10, 0.1, 1.0e-10, 0.0, 1.0])
    up_lims = np.array([1e3, 10, 1.0e-6, 2.0e-11, 10])
    init_guess = np.array([197.5, 1.87, 2.3e-8, 0.4e-11, 2])
    ydata = dark_df[['45CO2','47CO2','49CO2']].T.to_numpy().flatten()
    #print(ydata.shape, type(ydata))
    #to fit ydata, (3,N) array cannot be calculated, so fatten
    xdata = dark_df.index.to_numpy()#time array
    #print(xdata.shape, type(xdata))
    #j = f_dark(xdata, *init_guess)
    #print(j.shape, type(j))
    popt2, pcov2 = curve_fit(f_dark, xdata, ydata, \
                              p0 =init_guess,bounds = (low_lims, up_lims))
        
    dark_df.plot(y=['45CO2','47CO2','49CO2'], style ='o', markersize = 2)
    CO2s = f_dark(xdata,*popt2).reshape(3,-1)
    labels = ['C45', 'C47', 'C49']
    for i in range(3):
        plt.plot(xdata,CO2s[i,:], '--', label = labels[i]+'fit')
    show_plot()
    return popt2    
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

def calc_HCO3t(O2t, CO2t, HCO3_init = 0, O2_init = 0, CO2_init = 0, PQ = 1):
    """
    Estimate the [HCO3-] at a given time during light illumination
    Physiol. Plant, 90, 1994, Badger eat 1994, equation(8)
    this estimation is based on n(HCO3- extracellular) >> n(HCO3- cellular)
    Parameters
    ----------
    O2t : float or array
        [O2] at timepoint(s) t.
    CO2t : float or array
        [CO2] at timepoint(s) t.
    HCO3_init : float, optional
        Initial (dark before turning on light) HCO3- concentration.
    O2_init : float, optional
        Initial [O2]. The default is 0.
    CO2_init : float, optional
         Initial [CO2]. The default is 0.
    PQ : float, optional
        photosynthetic quotient. ∆O2 = ∆DIC * PQ. The default is 1.
        
    Returns
    -------
    HCO3t : [HCO3-] at given time t, as one point or an array

    """
    HCO3t = HCO3_init + CO2_init - (O2t - O2_init)/PQ - CO2t
    #here the PQ is assumed as 1, at the begining of a light period, the PQ is 
    #harder to estimate due to the relative higher DIC uptake and lower O2 
    #evolution. At steady state, the value will be closer to true value.
    return HCO3t

def calc_Uc(HCO3, CO2, dCO2dt):
    """
    calculate CO2 uptake to the cell
    Physiol. Plant, 90, 1994, Badger eat 1994, equation (4)
    similar to Hopkinson et al. 2011 PNAS, [S25] which has unit of mol/s
    here is to calculate the rate in µM/s(or similar), Uc_PNAS = Uc*Ve
    Parameters
    ----------
    HCO3 : float or 1-D array
        total [HCO3-]e.
    CO2 : float or 1-D array
        total [CO2]e.
    dCO2dt : float or 1-D array
        total d[CO2]/dt, MIMS signal.

    Returns
    -------
    Uc, CO2 uptake rate, unit µM/s.

    """
    global k1, k2
    Uc = k2 * HCO3 - k1 * CO2 - dCO2dt
    return Uc

def calc_Ub(dO2dt, Uc, PQ =1):
    """
    Calculate HCO3- upatke rate
    Physiol. Plant, 90, 1994, Badger eat 1994, equation (5)
    similar to Hopkinson et al. 2011 PNAS, [S26] which has unit of mol/s
    here is to calculate the rate in µM/s(or similar), Ub_PNAS = Ub*Ve
    Parameters
    ----------
    dO2dt : float or 1-D array
        O2 evolution rate, µM/s or similar unit.
        the proxy of CO2 fixation rate
    Uc : float or 1-D array
        CO2 uptake rate.

    Returns
    -------
    Ub, HCO3- upatke rate.

    """
    Ub = dO2dt/PQ -Uc#using PQ =1 in accordance with earlier work
    return Ub

def smt_curve(sig_df, window = 49, order = 1, PLOT = False):
    """
    smoothe a noisy signal curve using multi-step savgol_filter

    Parameters
    ----------
    sig_df : Pandas.Series
        the sigal to be smoothed
    window : int, optional
        maximum window size. The default is 49.
    order: int, optional
        savgol_filter order. The default is 1.
    PLOT: bool, optional
        whether or not to plot the smoothing process, the default is False

    Returns
    -------
    sig: a smoothed sigal as numpy.array

    """


    w = 3
    x= sig_df.index
    sig = sig_df.to_numpy()
    plt.plot(x, sig, label='raw signal')    
    while w<=window:
        sig = savgol_filter(sig, w, order)
        
        if PLOT:
            plt.plot(x, sig)
        w = w+2
    plt.plot(x, sig, label='smoothed sig. window='+str(window))    
    show_plot()
    return sig
 
def calc_do2dt_dco2dt(light_df, linear=False, \
                      window=33, order =1, C13_only = False):
    """
    

    Parameters
    ----------
    light_df : pd.DataFrame
        dark or light period df.
    linear : Bool, optional
        True or False.
    window : int, optional
        window size for salgol, odd number. The default is 33.
    order : int, optional
        salgol_filter parameter. The default is 1.
    C12: Bool, optional  

    Returns
    -------
    (dco2dt, do2dt, smt_CO2, smt_O2) in µM/s and µM.

    """
    global mims_do2dt
    if '47CO2' in light_df:
        CO2s = ['CO2','45CO2','47CO2','49CO2']
        if C13_only:
            CO2s = CO2s[1:]
        light_df['total_CO2'] = light_df[CO2s].sum(axis=1)
    else:
        CO2s = ['CO2','45CO2']
        if C13_only:
            CO2s = CO2s[1:]
        light_df['total_CO2'] = light_df[CO2s].sum(axis=1)
    #window = 33
    #order = 1
    if 'O2_Ar' in light_df:
        oxygen = 'O2_Ar'
    else:
        oxygen = 'O2'
    smt_CO2 = smt_curve(light_df['total_CO2'], window, order)
    smt_O2 = smt_curve(light_df[oxygen], window, order)
    fig1, ax1 = plt.subplots(sharex = True, dpi =300)
    light_df.plot(y=['total_CO2'],ax = ax1, legend = False, marker = 'o', \
                color = 'k', markersize = 2)
    ax1.plot(light_df.index, smt_CO2, 'g-')
    ax1.set_ylabel('$CO_2$/µM')
    ax1.legend(['total_$CO_2$','savgol_$CO_2$'], loc='upper left')
    ax2 = ax1.twinx()
    light_df.plot(y=[oxygen],ax = ax2,  legend = False, marker ='o',\
                color = 'b', markersize =2)
    ax2.plot(light_df.index, smt_O2,'m-')
    ax2.set_ylabel('$O_2$/µM')
    
    ax2.legend(['$O_2$','savgol_$O_2$'], loc='lower right')
    plt.show()
    plt.close()
    #the following estimate the first order derivative
    dco2 = savgol_filter(smt_CO2, window, order, deriv =1)
    
    dt = (light_df.index[-1]-light_df.index[0])/(light_df.index.size-1)
    dco2dt = dco2/dt
    if linear:#fit O2 evolution as a linear line        
        slope, intercept, r_value, p_value, std_err = \
            stats.linregress(light_df.index, light_df[oxygen].to_numpy())
        do2dt = np.zeros(light_df.index.size)#do2/dt
        do2dt = do2dt + slope -mims_do2dt #mims consumption correction
    else:
        do2 = savgol_filter(smt_O2, window, order, deriv =1)
        do2dt =do2/dt - mims_do2dt #mims consumption correction
    return (dco2dt, do2dt, smt_CO2, smt_O2, fig1)

def calc_PS_Urates(light_df,  HCO3_0, linear=True, window = 33, order =1,\
                   PQ_light = 1):
    """
    calculate the rate of photosynthesis as well as DIC uptake rates
    savgol_filter is used in the calculation
    Parameters
    ----------
    lihgt_df : pandas.DataFrame
        the df that holds O2, CO2 signal.
    PQ_light : float, optional
        the photosynthetic quotient of the light period

    Returns
    -------
    a DataFrame of following parameters
    smt_CO2:the smoothed total 13CO2 concentration;
    smt_HCO3: the smoothed estimation of [HCO3-];
    smt_O2:   the smoothed O2 concentration
    PS_rate:  the photosynthetic rate as O2 evolution rate per cell.
    Uc:       the CO2 uptake rate, mol/s/cell or similar units
    Ub:       the HCO3- uptake rate, mol/s/cell or similar units
    srf_CO2:
    srf_HCO3:
    Uc_M:
    Ub_M:
    """
    global Ve, Vs, N, fc_BL, fb_BL, kSF, kSR
    dco2dt, do2dt, smt_CO2, smt_O2, fig1 = \
        calc_do2dt_dco2dt(light_df, linear, window, order)
    """
    note that savgol returns the derivative assuming spacing as indexing
    which means dt = 1, since dt in the experiment should be:
    dt = (light_df.index[-1]-light_df.index[0])/(light_df.index.size-1)
    or 3.58 seconds, this dt may vary if experiments are carried differently
    
    so real do2dt needs corretion from dt spacing, 
    """
    fig2, ax2 = plt.subplots(dpi = 300)
    ax2.plot(light_df.index, dco2dt,'go', label = 'd$[CO_2]$/dt')
    ax2.set_ylabel('d$[CO_2]$/dt', color = 'g')
    ax2.legend(loc='upper left')
    ax3 = ax2.twinx()
    ax3.plot(light_df.index, do2dt, 'mo', label = 'd$[O_2]$/dt')
    ax3.set_ylabel('d$[O_2]$/dt', color = 'm')
    ax3.legend(loc='lower right')
    plt.show()
    plt.close()
    
    
    
    smt_HCO3 = calc_HCO3t(smt_O2, smt_CO2, HCO3_init = HCO3_0, \
                     O2_init=smt_O2[0], CO2_init = smt_CO2[0], PQ = PQ_light)
    Uc = calc_Uc(smt_HCO3, smt_CO2, dco2dt)
    Ub = calc_Ub(do2dt, Uc, PQ = PQ_light)
    fig3 = plt.figure(dpi = 300)
    plt.plot(light_df.index, do2dt,'m-')
    plt.plot(light_df.index, Uc, 'g-')
    plt.plot(light_df.index, Ub, 'r-')
    plt.plot(light_df.index, np.zeros(light_df.index.size), 'k--')
    plt.legend(['$O_2$ evolution', '$CO_2$ uptake', '$HCO_3^-$ uptake'])
    plt.title('Phothosynthesis rate and DIC uptake')
    plt.xlabel('time/s')
    plt.ylabel('rate: µM/s')
    plt.show()
    plt.close()
    
    #now calculate surface values
    srf_CO2 = smt_CO2 - Uc*Ve/N/fc_BL
    srf_HCO3= smt_HCO3- Ub*Ve/N/fb_BL
    dco2_srf = savgol_filter(srf_CO2, window, order, deriv =1)
    dt = (light_df.index[-1]-light_df.index[0])/(light_df.index.size-1)
    dco2dt_srf = dco2_srf/dt

    Uc_M = Uc * Ve/N/Vs + kSR*srf_HCO3 - kSF*srf_CO2 - dco2dt_srf
    Ub_M = calc_Ub(do2dt*Ve/N/Vs, Uc_M, PQ = PQ_light)
    #note that the Ve is for the whole system, so the impact of Uc(bulk) and
    #do2dt should be normalized to each cell, divided by N
    #so far the Uc, Ub are calculated for N cells and Uc_M, Ub_M for single cell.
    #that is, Uc, Ub in unit of µM/(N cells)/s, while Uc_M, Ub_M is in µM/cell/s
    
    #Later these values are all coverted to mol/cell/s
    
    # figx = plt.figure(dpi = 300)
    # plt.plot(light_df.index, do2dt*Ve/N/Vs,'m-')
    # plt.plot(light_df.index, Uc_M, 'g-')
    # plt.plot(light_df.index, Ub_M, 'r-')
    # plt.plot(light_df.index, np.zeros(light_df.index.size), 'k--')
    # plt.legend(['$O_2$ evolution', '$CO_2$ uptake', '$HCO_3^-$ uptake'])
    # plt.title('Phothosynthesis rate and Surface DIC uptake')
    # plt.xlabel('time/s')
    # plt.ylabel('rate: µM/s')
    # plt.show()
    # plt.close()
    
    Uc_M = Uc_M*Vs*1e-9 #covert to mol/cell/s
    Ub_M = Ub_M*Vs*1e-9 #covert to mol/cell/s
    
    
    Uc = Uc*Ve*1e-9/N #covert to mol/cell/s
    Ub = Ub*Ve*1e-9/N #covert to mol/cell/s
    do2dt = do2dt*Ve*1e-9/N #covert to mol/cell/s
    fig4 = plt.figure(dpi = 300)
    plt.plot(light_df.index, do2dt,'m-')
    plt.plot(light_df.index, Uc, 'g-')
    plt.plot(light_df.index, Ub, 'r-')
    plt.plot(light_df.index, np.zeros(light_df.index.size), 'k--')
    plt.legend(['$O_2$ evolution', '$CO_2$ uptake', '$HCO_3^-$ uptake'])
    plt.title('Phothosynthesis rate and DIC uptake')
    plt.xlabel('time/s')
    plt.ylabel('rate: mol/s/cell')
    plt.show()
    plt.close()
    
    result_df = pd.DataFrame({'smt_CO2':smt_CO2,\
                              'smt_HCO3':smt_HCO3,\
                              'smt_O2': smt_O2,\
                              'PS_rate': do2dt,\
                                  'Ub': Ub,\
                                  'Uc': Uc,\
                              'srf_CO2':srf_CO2,\
                              'srf_HCO3':srf_HCO3,\
                                'Uc_M':Uc_M,\
                                'Ub_M':Ub_M}, index = light_df.index)

    ax=result_df.plot(y=['Ub_M','Ub','PS_rate','Uc','Uc_M'])
    ax.plot(light_df.index, np.zeros(light_df.index.size), 'k--')
    ax.set_ylabel('rate: mol/cell/s')
    plt.title('PS rate and DIC use bulk vs surface')
    fig5 = ax.get_figure()  
    show_plot()
    
    figs = [fig1, fig2, fig3, fig4, fig5]
    return (result_df, figs)

def cal_Kks(pH, Tc=25, S=35):
    """
    Calculate K1, K2, k1, k2

    Parameters
    ----------
    pH : scalar or array
        DESCRIPTION.
    Tc : scalar, optional
        temperature in C. The default is 25.
    S : scalar, optional
        salinity in per mil. The default is 35.

    Returns
    -------
    K1, K2, k1, k2 with a given pH

    """
    T = Tc + 273.15
    H = 10**(-pH)
    pK1 = 3633.86/T - 61.2172 + 9.6777*np.log(T) - 0.011555*S + 1.152e-4 *S**2
    K1 = 10**-pK1          
    #K1 equilibrium constant between CO2 and HCO3 Lueker, Dickson, Keeling Mar Chem. 2000.
    pK2 = 471.78/T + 25.929 - 3.16967 * np.log(T) - 0.01781*S + 1.122e-4 *S**2
    K2 = 10**-pK2
    #K2 equilbrium constant from Lueker, Dickson, Keeling Mar Chem 2000
    Kw = np.exp(148.96502 - (13847.26 / T) - 23.6521 * np.log(T)\
            + (S**0.5)*((118.67 / T) - 5.977 + 1.0495 * np.log(T)) - 0.01615 * S) 
    bfrac_e = (1/ (1+ K2/H))       # fraction of "B" pool that is HCO3- in extracellular solution
    kp1 = np.exp(1246.98 - 6.19E4 / T - 183 * np.log(T))  # CO2 + H2O -> H+ + HCO3- Johnson 1982 as presented in Zeebe and Wolf-Gladrow
    kp4 = 4.70E7 * np.exp(-23.2 / (8.314E-3 * T))  # CO2 + OH- -> HCO3- Johnson 1982
    # km1 = kp1/K1            # H+ + HCO3- -> CO2 + H2O
    # km4 = kp4*Kw/K1         # HCO3- -> CO2  + OH-
    kuf = kp1 + kp4 * (Kw/H)      #CO2 hydration rate in bulk solution
    kur = bfrac_e * kuf * (H/K1)  #HCO3- dehyration rate in bulk solution
    K_dict = {'pK1':pK1,'pK2':pK2, 'K1': K1, 'K2': K2, 'k1': kuf, 'k2': kur}        
    return K_dict 

def calc_RS(dark_df, linear=True, window = 33, order =1):
    """
    Calculate Respiration Rate after light off

    Parameters
    ----------
    dark_df : pd.DataFrame
        dark data.
    linear : BOOL, optional
        DESCRIPTION. The default is True.
    window : int, optional
        Savgol_filter window. The default is 33.
    order : int, optional
        Savgol_filter order. The default is 1.

    Returns
    -------
    result_df : TYPE
        DESCRIPTION.
    fig : matplot object
        Figure.

    """
    dco2dt, do2dt, smt_CO2, smt_O2, fig1 = \
        calc_do2dt_dco2dt(dark_df, linear, window, order)

    do2dt = do2dt*Ve*1e-9/N    
    result_df = pd.DataFrame({'RS_rate': do2dt}, index = dark_df.index)
    fig = fig1
    return (result_df, fig)

def calc_HCO3_dark(light_df, linear=True, window = 33, order =1,\
                   C13only = False, fO2_CO2 = 1, eq = False):
    """
    Calculate HCO3- concentrations in the dark

    Parameters
    ----------
    light_df : dataframe 
        before a light period.
    linear : BOOL, optional
        Linear fit for O2. The default is True.
    window : int, optional
        Savgol_filter window. The default is 33.
    order : int, optional
        Savgol_filter order. The default is 1.
    C13only: BOOL, optional
        whether or not only considering 13C
    fO2_CO2: numeric, optional
        fraction of O2 consumption rate to CO2 increase rate in the bulk
        there could be a CO2 --> HCO3- conversion at the cell surface before
        the realse of DIC into the bulk in the dark. 
        The default is 1, as described in Badger 1994
    eq: BOOL, optional
        whether or not CO2 reaches equilibrium with HCO3- at cell surface
        if True, this will change fO2_CO2 to k2/(k1+k2),
        otherwise, it does not impact fO2_CO2, as The default is False

    Returns
    -------
    HCO3_ave, averaged HCO3 concentration

    """
    global k1, k2
    dco2dt, do2dt, smt_CO2, smt_O2, fig1 = \
        calc_do2dt_dco2dt(light_df, linear, window, order, C13only)
    fig2, ax2 = plt.subplots(dpi = 300)
    ax2.plot(light_df.index, dco2dt,'go', label = 'd$[CO_2]$/dt')
    ax2.set_ylabel('d$[CO_2]$/dt (µM/s)', color = 'g')
    ax2.legend(loc='upper left')
    ax3 = ax2.twinx()
    ax3.plot(light_df.index, do2dt, 'mo', label = 'd$[O_2]$/dt')
    ax3.set_ylabel('d$[O_2]$/dt (µM/s)', color = 'm')
    ax3.set_xlabel('time/s')
    ax3.legend(loc='lower right')
    plt.show()
    plt.close()
    
    
    #the following eq7 is from Badger 1994
    #HCO3_dark = (dco2dt + k1*smt_CO2 + do2dt)/k2
    #without eCA, the above equation may be true
    #with eCA, one may add a factor of fO2_CO2 to see the impact of 
    #surface equilibrium on HCO3- calculation:
    #HCO3_dark = (dco2dt + k1*smt_CO2 + do2dt*fO2_CO2)/k2
    if eq:
        fO2_CO2 = k2/(k1+k2)    
    if C13only:
        HCO3_dark = (dco2dt + k1*smt_CO2 + do2dt*0.01109 *fO2_CO2)/k2
        #0.01109 is the nature abundance of 13C
    else:
        HCO3_dark = (dco2dt + k1*smt_CO2 + do2dt *fO2_CO2)/k2
    
    #note that resp rate = -do2dt, and leaking rate is not considered
    #so this equation does not apply right after light is turned off
    fig3, ax3 = plt.subplots(dpi = 300)
    ax3.plot(light_df.index, HCO3_dark,'mo', label = '$[HCO_3^-]$')
    ax3.set_ylabel('$[HCO_3^-] (µM)$', color = 'm')
    if C13only:
        ax3.set_ylabel('$[H^{13}CO_3^-] (µM)$', color = 'm')
    ax3.legend(loc='lower right')
    plt.show()
    plt.close()
    
    HCO3_ave = HCO3_dark[-10:].mean()
    return HCO3_ave

def calc_DIC_dark(light_df, linear=True, window = 33, order =1,\
                  C13only = False, fO2_CO2 = 1, eq = False):
    """
    Calculate CO2 and HCO3- concentrations in the dark

    Parameters
    ----------
    light_df : dataframe 
        before a light period.
    linear : BOOL, optional
        Linear fit for O2. The default is True.
    window : int, optional
        Savgol_filter window. The default is 33.
    order : int, optional
        Savgol_filter order. The default is 1.
    C13only: BOOL, optional
        whether or not only considering 13C
    fO2_CO2: numeric, optional
        fraction of O2 consumption rate to CO2 increase rate in the bulk
        there could be a CO2 --> HCO3- conversion at the cell surface before
        the realse of DIC into the bulk in the dark. 
        The default is 1, as described in Badger 1994
    eq: BOOL, optional
        whether or not CO2 reaches equilibrium with HCO3- at cell surface
        if True, this will change fO2_CO2 to k2/(k1+k2),
        otherwise, it does not impact fO2_CO2, as The default is False

    Returns
    -------
    CO2, last 10 smt CO2 values' average
    HCO3_dark, last 10 smt HCO3 values' average

    """
    global k1, k2
    dco2dt, do2dt, smt_CO2, smt_O2, fig1 = \
        calc_do2dt_dco2dt(light_df, linear, window, order, C13only)
    fig2, ax2 = plt.subplots(dpi = 300)
    ax2.plot(light_df.index, dco2dt,'go', label = 'd$[CO_2]$/dt')
    ax2.set_ylabel('d$[CO_2]$/dt (µM/s)', color = 'g')
    ax2.legend(loc='upper left')
    ax3 = ax2.twinx()
    ax3.plot(light_df.index, do2dt, 'mo', label = 'd$[O_2]$/dt')
    ax3.set_ylabel('d$[O_2]$/dt (µM/s)', color = 'm')
    ax3.set_xlabel('time/s')
    ax3.legend(loc='lower right')
    plt.show()
    plt.close()
    
    
    #the following eq7 is from Badger 1994
    #HCO3_dark = (dco2dt + k1*smt_CO2 + do2dt)/k2
    #without eCA, the above equation may be true
    #with eCA, one may add a factor of fO2_CO2 to see the impact of 
    #surface equilibrium on HCO3- calculation:
    #HCO3_dark = (dco2dt + k1*smt_CO2 + do2dt*fO2_CO2)/k2
    if eq:
        fO2_CO2 = k2/(k1+k2) 
    if C13only:
        HCO3_dark = (dco2dt + k1*smt_CO2 + do2dt*0.01109 *fO2_CO2)/k2
        #0.01109 is the nature abundance of 13C
    else:
        HCO3_dark = (dco2dt + k1*smt_CO2 + do2dt *fO2_CO2)/k2
    
    #note that resp rate = -do2dt, and leaking rate is not considered
    #so this equation does not apply right after light is turned off
    fig3, ax3 = plt.subplots(dpi = 300)
    ax3.plot(light_df.index, HCO3_dark,'mo', label = '$[HCO_3^-]$')
    ax3.set_ylabel('$[HCO_3^-] (µM)$', color = 'm')
    if C13only:
        ax3.set_ylabel('$[H^{13}CO_3^-] (µM)$', color = 'm')
    ax3.legend(loc='lower right')
    plt.show()
    plt.close()
    
    HCO3_ = HCO3_dark[-10:].mean()
    CO2 = smt_CO2[-10:].mean()
    return (CO2, HCO3_)

def MIMS_resp(MIMS_df):
    if 'O2_Ar' in MIMS_df:
        oxygen = 'O2_Ar'
    else:
        oxygen = 'O2'
    MIMS_do2dt, intercept, r_value, p_value, std_err = \
                stats.linregress(MIMS_df.index, MIMS_df[oxygen].to_numpy())
    t_MIMS = MIMS_df.index.to_numpy()
    plt.plot(t_MIMS, MIMS_df['O2_Ar'].to_numpy(), 'bo', label = 'data')
    plt.plot(t_MIMS, MIMS_do2dt*t_MIMS+intercept, 'g--', label = 'fit')
    plt.ylabel('$O_2$/µM')
    plt.xlabel('time/s')
    show_plot()
    return MIMS_do2dt
#### #end of functions, start of analysis code#######

G = np.zeros((4,3))#matrix S52 
#Hopkinson et al. www.pnas.org/cgi/content/short/1018062108
for i in range(3):
    G[i,i] = 1
H = np.zeros((3,4))# matrix S54 
for i in range(3):
    H[i, i] = (3-i)/3
    H[i, i+1] = (i+1)/3

#### load data: calculated concentrations####
f_name = 'Fcyl_3C_S40__exp02_c.csv'
data = pd.read_csv(f_name, index_col=0)
data.index.rename('time/s', inplace=True)


#use 49CO2 to determine time of addition of NaHCO3
CO2_49 = data['49CO2'].iloc[:1200].to_numpy()
dC49 = (CO2_49[10:]+1)/(CO2_49[:-10]+1)
# plt.plot(data.index[:-10],dC49)
t0 = np.argwhere(dC49==dC49.max())[0][0]
t1 =np.argwhere(dC49[t0+60:]==dC49[t0+60:].min())[0][0]+t0+60
plot_df(data.iloc[t0:t1])
plot_df(data.iloc[t0-15:t0+10])
# t0 = t0+9
if 'O2_Ar' in data:#correction using Ar signal    
    for col in data.columns[-4:-1]:
        data[col] = data[col] - data[col].iloc[t0-60:t0].mean()
else:
    for col in data.columns[-3:]:
        data[col] = data[col] - data[col].iloc[t0-60:t0].mean()
# t1 = t1-5
passive = data.iloc[t0:t1]
passive.plot(y=['45CO2','47CO2','49CO2'], style ='o')
plot_df(passive.iloc[-30:])

t1 = t1-2
passive = data.iloc[t0:t1]
passive.plot(y=['45CO2','47CO2','49CO2'], style ='o')

#### HCO3- hydration simulation with approximate k1, k2 values and fit
# passive.iloc[-20:].plot(y=['45CO2','47CO2','49CO2'], style ='o')

#The fit_passive function solve the apparent values for k1 and k2
#The F18O is critical, should we solve this or measure F18O
k1,K_bulk, DIC0, f180 = fit_passive_3(passive.iloc[:])
print(k1,K_bulk, DIC0, f180)
#K_bulk is the uncatalyzed equilibrium constant.
k2 = k1/K_bulk
# k1, k2 = fit_passive(passive, F18O = 0.90, F13C = 0.985) 
# print(k1,k2)
#if 18Oand 13C decreases overtime, adjust those values for better fit
#k1, k2 are kuf, kur calculated from fitting data

#### Add cells, fc, fb, kcf... ####
#the following calculations involve cells; some constants are assigned first
N = 4.40e8*0.020#total number of cells
Ve = 1225*1e-3#rxn cuvette volume in mL
Vc = 41.6e-12#cellular volume in mL, measured average volume at 3°C
#measured average effective diameter ~4.30µm
Vs = 4/3 * np.pi * (2.25**3-2.15**3) * 1e-12
#assumming 0.1 µm surface space, unit in mL
r_cell = 2.15e-4 #unit in cm
fc_BL, fb_BL = 4*np.pi * r_cell * Di(3, 40)[:2]
#Di(3, 40) gives an array of CO2, HCO3- CO3-- diffusion coefficient
#at 3 °C, S=40, fc_BL, fb_BL in cm3/s, ref: Hopkinson 2013 et al Plant Physiol
fc, fb = 1.51e-9, 1.38e-12 #these values are calculated from exp. with AZ
fc_M, fb_M = 1/(np.array([1/fc, 1/fb]) - np.array([1/fc_BL, 1/fb_BL]))

# kcf, kcr = 194.2, 2.04

#carving the time points out for events after t1
"""The t0, t1, t2, t3, t4, t5, (and more if experiments are different)
    are the timepoint as index value of the following events:
    13C-18O-NaHCO3 addition, cell addition, light on, off, on, off,...
    If the approximate time difference between two time points is known,
    one can change the np.argwhere to specify the time window for event
    time point identification, for instance:
    t[i+1] = np.argwhere(dC45 == dC45[ti: ti+within_known_points].max())[0][0]
"""
#dC45 = conv_a_sig(data['45CO2'], diff=np.array([1]*10+[0]+[-1]*10))
dC47 = conv_a_sig(data['49CO2'], size=63, diff=np.array([1]*30+[0]+[-1]*30))
#the conv_a_sig will plot the result
#I found dC47 is useful in this experiment,
#one may unquote the dC45, dC49 line to see which is better for identify t
#dC49 = conv_a_sig(data['49CO2'], diff=np.array([1]*10+[0]+[-1]*10))
t = [t0, t1]# store time points of event in t list
ti=t1
for i in range(4):
    dC47[:ti+80] =0
    if i%2 == 0:
        ti = np.argwhere(dC47==dC47[:335+ti].min())[0][0]
    else:
        ti = np.argwhere(dC47==dC47[:700+ti].max())[0][0]
        #700+ti is to find max within 20 min(light on to light off),
        #16.77 points per min
    t.append(ti)    
    plt.plot(dC47, label = 'after t'+str(i+1))
    plt.title('convovled signal of 49CO2')    
show_plot()
# for i in range(len(t)-1):
#     plot_df(data.iloc[t[i]:t[i+1]])
# uncomment to see how the t list slices the data


#### Calculate fc, fb, kcf... ####
sol = sim_passive([data.index[t0],data.index[t1]], \
                  k1, k2, DIC_t = DIC0, f18O = f180, f13C=0.985)
y01 = sol.y[:,-1]*1.205/Ve ##1.205 is the volume before adding cells
y0s = y01#assuming fast exchange between bulk and surface, 
#otherwise the ksf and ksr cannot be used for differential equations
y0in = np.zeros(7)#intracelluar
#y0in[:3] = y01[:3]#assumming fast CO2 diffusion to cells
y01 = np.concatenate((y01, y0in, y0in))


dark0 = data.iloc[t[1]:t[2]-20]#there may be a delay due mixing
plot_df(dark0, False, False)
plot_totalC(dark0)

kSF, kcf, Kcyto = fit_dark_SF(dark0, y01)
kSR = kSF*k2/k1
kcr = kcf/Kcyto
print('kSF, kcf, kcr, Kcyto= ', kSF, kcf, kcr, Kcyto)


ctot = data[['CO2','45CO2','47CO2','49CO2']].sum(axis = 1)
ctot.plot()


#### first light period
light1 = data.iloc[t[2]-20:t[3]-12]
plot_df(light1, True, True)

#dark period following light 1
dark1 = data.iloc[t[3]-12:t[4]-15]
plot_df(dark1, True, True)

#### calculate MIMS O2 consumption ####

#even though O2 is Ar corrected, the consumption of O2 seems faster than Ar
#A MIMS_do2dt can be calculated for real rate adjustment, assuming constant
mims_df = data.iloc[t0+17:t1-17]#avoiding two ends of the NaHCO3 hydration part
plot_df(mims_df)

mims_do2dt = MIMS_resp(mims_df)#unit in µM/s

#### calculate the light1 rates

CO2_dark0, HCO3_dark0 = calc_DIC_dark(dark0.iloc[85:],linear= True, window=41, eq = True)
#two array of values

DIC_dark0 = CO2_dark0 + HCO3_dark0
l1_df, Figs_l1 = calc_PS_Urates(light1.iloc[:-17], HCO3_dark0, False, window = 51)

l1_df['cyt_CO2'] = l1_df['srf_CO2'] - l1_df['Uc_M']/fc_M*1e9
l1_df.to_csv(f_name[:-14]+'_1st_light_noAZ.csv')


#### prepare for light 2 calculations ####    
t4_CO2 = ctot.iloc[t[4]-25:t[4]-15].mean()
t4_O2 = data['O2_Ar'].iloc[t[4]-25:t[4]-15].mean()
t2_O2 = data['O2_Ar'].iloc[t[2]-30:t[2]-20].mean()
delta_O2 = t4_O2 - t2_O2 - (data.index[t[4]-15] -data.index[t[2]-20])*mims_do2dt
t4_HCO3 = HCO3_dark0 - delta_O2/1 + CO2_dark0 - t4_CO2

DIC_t4 = t4_HCO3 + t4_CO2
#∆HCO3 = -∆O2/PQ - ∆CO2, PQ is the photosynthetic quotian, here 1 is used

# t4_HCO3 = calc_HCO3_dark(dark1[100:], linear = False, window =33)

#t4_HCO3 = k1*ctot.iloc[t[4]].mean()/k2
light2 = data.iloc[t[4]-15:t[5]-4]
plot_df(light2, True, True)
plot_totalC(light2)

dark2 = data.iloc[t[5]-4:]
plot_df(dark2)

l2_df, Figs = calc_PS_Urates(light2.iloc[:-17], t4_HCO3, False, window=51)

l2_df['cyt_CO2'] = l2_df['srf_CO2'] - l2_df['Uc_M']/fc_M*1e9
l2_df.to_csv(f_name[:-14]+'_2nd_light_noAZ.csv')



#### calculate steady state rates ####
l1 = l1_df.iloc[-85:].mean()
l1['smt_DIC'] = l1.smt_CO2 + l1.smt_HCO3
l1['DICo'] = DIC_dark0

d1, rsfig1 = calc_RS(dark1)
l1['RS_rate'] = d1['RS_rate'].mean()
Figs_l1.append(rsfig1)

l2 = l2_df.iloc[-85:].mean()
l2['smt_DIC'] = l2.smt_CO2 + l2.smt_HCO3
l2['DICo'] = DIC_t4

d2, rsfig2=calc_RS(dark2)
l2['RS_rate'] = d2['RS_rate'].mean()
Figs.append(rsfig2)

pdf_name = f_name[:-14]+'_plots_noAZ.pdf'
pp = PdfPages(pdf_name)
for afig in Figs_l1+Figs:
    pp.savefig(afig)
pp.close()

l1['light_on'] = '1st'
l2['light_on'] = '2nd'

cols = ['Date','S','T/°C','k1','k2','K_blk','ksf','ksr', 'kcf', 'kcr','K_cyt',\
        'fc','fb','fc_BL','fb_BL','fc_M','fb_M','N_total','Vs','Vc','Ve',\
            'AZ','MIMS_do2dt']
paras = [f_name.split()[0][-8:], 40, 3, k1, k2, K_bulk, kSF, kSR, kcf, kcr, Kcyto,\
         fc, fb, fc_BL, fb_BL, fc_M, fb_M, N, Vs, Vc, Ve,\
             'without', mims_do2dt*Ve/N*1e-9]

para_s = pd.Series(data = paras, index=cols)
l1_s = pd.concat([para_s, l1])
l2_s = pd.concat([para_s, l2])
smr_noAZ = pd.DataFrame(data=[l1_s,l2_s], columns = l1_s.index)
