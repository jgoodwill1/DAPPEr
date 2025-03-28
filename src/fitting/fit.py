"""
Created on March 25, 2023

@author: Josh Goodwill
"""

import numpy as np
from scipy.optimize import curve_fit

'''
def gen_fit(x, y):
Generates fit of linear and sqrt portion of IV arrays. Refer to
data_processing to identify filtering of IV arrays. scipy.optimize.curve_fit
used to fit non-linear least squares to fit


input:
V_arr (array); Voltage array
I_arr (array); Current array
proc* (boolean); Whether to process the arrays for better fitting

output:
V_fit (np.array); linear space between min and max voltages from V_arr
model(t, *popt) (np.array); fitted model of I_arr using Voltage range and 
                            pointer to popt
popt(array); optimal values for parameters
V0   (float) = popt[0]; floating potential [V]
VP   (float) = popt[1]; plasma potential [V]
m1   (float) = popt[2]; slope of linear fit
b    (float) = popt[3]; y-intercept of linear fit
ne   (float) = popt[4]; electron density [cm^{-3}]
Te   (float) = popt[5]; electron temperature [K]


Intial guesses/[bounds]:
    V0    = 2.0         ;       [-3         :           3]
    VP    = 2.0         ;       [2          :           5]
    m1    = 10          ;       [-1000      :        1000]
    y_int = 80          ;       [-1000      :        1000]
    ne    = 5*10^10     ;       [5          :  5*(10**15)]
    Te    = 1500        ;       [300        :        5000]

pcov(2D np.array); covariance of popt array
'''

def gen_fit(V_arr, I_arr, proc = False):
    V_proc = []
    I_proc = []
    if proc == True:
        V_proc, I_proc = data_processing(V_arr, I_arr)
    else:
        V_proc, I_proc = V_arr, I_arr
    g = [ 2,   2,    10,    80, 5*(10**10),  1500]    #intial guess
    b = [
        (-3,   2, -1000, -1000,  5*(10**1),   300),
        ( 3,   5,  1000,  1000, 5*(10**15),  5000)
        ] #bounds
    popt, pcov = curve_fit(model,V_proc, I_proc, g, bounds=b)
    V_fit = np.linspace(min(V_proc),max(V_proc), num = 300) #Voltage array processed for fit
    return V_fit, model(V_fit,*popt), popt, pcov

'''
def model(V_proc, VP, m1, y_int, ne, Te, V0):
Model to fit ion, *transition*, electron saturation regions

input:
V_proc (np.array); processed V_arr for fitting
V0 
VP (float); plasma potential [V]
m1 (float); slope of linear fit
b  (float); y-intercept of linear fit
ne (float); electron density [cm^{-3}]
Te(float); electron temperature [K]

'''

def model(V, V0, VP,  y_int, m1, ne, Te):
    I = np.zeros(len(V))
    I[V <= VP] = lin(V[V <= VP], m1, y_int) - lin(VP, m1, y_int)
    # Vtrans = (V > VP)
    # I[Vtrans] = exp(V[Vtrans], ne, te, VP) + b
    I[V >  VP] = sqrt(V[V >  VP], ne, Te, V0) - sqrt(VP, ne, Te, V0) + y_int
    return I

def lin(x, m, y_int): #linear--full model square root
    return m * x + y_int

# def exp_fit(x, a, etemp, Vf): #exponential fit
#     q_e = 1.602 * 10**-19 #electron charge [C]
#     K_b = 1.381 * 10**-23 #boltzmann constant [m^2*kg/(s^2*K)]   
#     k = q_e / (K_b * Te)
#     return a * np.exp(k * (x - Vf))

'''
def sqrt_fit(x, ne, Te, V0)
Uses Eq. 1.4 and 1.5 of https://digitalcommons.usu.edu/etd/274

input:
x (np.array); V_proc [V]
ne (float); electron density [cm^{-3}]
Te(float) ; electron temperature [K]
V0 (float); plasma potential [V]

output:
I (float); current fit [nA]
'''
def sqrt(x, ne, Te, V0):# square root
    q_e = 1.602 * 10**-19 #electron charge [C]
    K_b = 1.381 * 10**-23 #boltzmann constant [m^2*kg/(s^2*K)]   
    m_e = 9.109 * 10**-31 #electron mass [kg]              
    R = (3./16.) * 0.0254 #radius of probe [cm?]
    L = (3.25) * 0.0254 #length of probe [cm?]
    A = 2. * np.pi * R * L + np.pi * (R ** 2) #top and length area of cylinder [cm^2]

    k = q_e / (K_b * Te)
    I0 = ne * q_e * np.sqrt(K_b * Te / (2. * np.pi * m_e)) * A / (10**-9)
    return I0 * np.sqrt(1. + k*(x - V0))

'''
def data_processing(V, I)
Processes Voltage and Current arrays for fitting algorithm.
Removes nan values and sets them to 0.

input:
V (np.array); Normal Voltage array [V]
I (np.array); Normal Current array [nA]

output:
V_proc (np.array); processed Voltage array
I_proc (np.array); processed Current array
'''

def data_processing(V, I):#remove data points below -2 and above the peak to reduce datapoints going to the fitting routine
    V_proc = np.nan_to_num(x_raw, nan=0.0)
    I_proc = np.nan_to_num(y_raw, nan=0.0)
    return V_proc, I_proc