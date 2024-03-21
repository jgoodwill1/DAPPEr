"""
Created on July 31 12:53:14 2021

@author: Josh Goodwill
"""

import numpy as np
from scipy.optimize import curve_fit

'''
def generate_fit(x, y):
Generates fit of linear and sqrt portion of IV arrays. Refer to
data_processing to idenify filtering of IV arrays. scipy.optimize.curve_fit
used to fit non-linear least squares to fit

input:
V_arr (array); Voltage array
I_arr (array); Current array

output:
V_fit (np.array); linear space between min and max voltages from V_arr
model(t, *popt) (np.array); fitted model of I_arr using Voltage range and 
                            pointer to popt
popt(array); optimal values for parameters
  V_f = popt[0]; Plasma potential
  I_p = popt[1]; Current at Floating 
  
pcov(2D np.array); covariance of popt array
'''

def gen_fit(V_arr, I_arr):
    # V_proc, I_proc = data_processing(V_arr, I_arr)
    V_proc, I_proc = (V_arr[V_arr < 4], I_arr[V_arr < 4])
    guess = [0.6,-14,80, 5*(10**10),1000,-0.5]    #intial guess
    b = ((-3,-np.inf,-np.inf,0,0,-3),(3,np.inf,np.inf,np.inf,10000,3)) #bounds
    popt, pcov = curve_fit(model, V_proc, I_proc, guess, bounds = b)
    V_fit = np.linspace(min(V_proc),max(V_proc), num = 50) #Voltage array processed for fit
    return V_fit, model(V_fit,*popt), popt, pcov, V_proc, I_proc

'''
def model(V_proc, VP, m1, b, ne, etemp, V0):
Model to fit ion, *transition*, electron saturation regions

input:
V_proc (np.array); processed V_arr for fitting
Vf (float); floating potential
m1 (float); slope of linear fit
b  (float); y-intercept of linear fit
ne (float); electron density [cm^{-2}]
etemp(float); electron temperature [K]
VP (float); negative plasma potential [V]
'''

def model(V_proc, Vf, m1, b, ne, etemp, VP):
    I_fit = np.zeros(len(V_proc))
    #Linear fit for ion saturation
    I_fit[V_proc <= Vf] = lin_fit(V_proc[V_proc <= Vf], m1, b) 
    #Sqrt fit for electron saturation
    I_fit[V_proc > Vf] = sqrt_fit(V_proc[V_proc > Vf], ne, etemp ,VP) \
                         - sqrt_fit(Vf, ne , etemp, VP)
    return I_fit

def lin_fit(x, m, y_int): #linear--full model square root
    return m * x + y_int

'''
def sqrt_fit(x, ne, etemp, V0)
Uses Eq. 1.4 and 1.5 of https://digitalcommons.usu.edu/etd/274

input:
x (np.array); V_proc [V]
ne (float); electron density [cm?]
etemp(float); electron temperature [K]
VP (float); negative plasma potential [V]

output:
I (float); current fit [nA]
'''
def sqrt_fit(x, ne, etemp, VP):# square root
    q_e = 1.602 * 10**-19 #electron charge [C]
    K_b = 1.381 * 10**-23 #boltzman constant [m^2*kg/(s^2*K)]   
    m_e = 9.109 * 10**-31 #electron mass [kg]              
    R = (3./16.) * 0.0254 #radius of probe [cm?]
    L = (3.25) * 0.0254 #length of probe [cm?]
    A = 2. * np.pi * R * L + np.pi * (R ** 2) #top and length area of cylinder [cm^2]

    k = q_e / (K_b * etemp)
    I0 = ne * q_e * np.sqrt(K_b * etemp / (2. * np.pi * m_e)) * A / (10**-9)
    return I0 * np.sqrt(1. + k*(x + VP))

'''
def data_processing()
'''

def data_processing(x_raw,y_raw):#remove data points below -2 and above the peak to reduce datapoints going to the fitting routine
    max_ind = np.argmax(y_raw) #find index of maximum y
    if (max_ind.size != 1):
        min_max_val = np.argmin(x_raw[max_ind]) #find minimum x of max y
    else:
        min_max_val = max_ind
    ind_high = np.where(x_raw > x_raw[min_max_val] + 0.1) # find indexs of x greater than 0.1 more than max y
    ind_low = np.where(x_raw < 0) #find indexs of x less than 0
    ind_rem = np.concatenate((ind_high,ind_low),axis=None)
    x = np.delete(x_raw,ind_rem) #remove from x
    y = np.delete(y_raw,ind_rem) #remove from y
    return x,y

def remove_outliers(time, temp, density): #finds and removes outliers
    ind_low_t = np.where(np.array(temp) < 0) #finds negative temperatures
    ind_high_t = np.where(np.array(temp) > 3 * np.percentile(temp,95)) #finds temperatures more than triple the 95 percentile
    ind_low_d = np.where(np.array(density) < 0) #finds negative density
    ind_high_d = np.where(np.array(density) > 3 * np.percentile(density,95)) #finds densities more than triple the 95 percentile
    ind_rem = np.unique(np.concatenate((ind_high_t,ind_low_t,ind_low_d,ind_high_d),axis=None))
    time_n = np.delete(time,ind_rem)
    temp_n = np.delete(temp,ind_rem)
    density_n = np.delete(density,ind_rem)
    return time_n,temp_n,density_n

def findIndexs(mypackets): #find indexs of various packet types
    sensor = []
    medium = []
    large = []
    burst = []
    i = 0
    #pcktType
    typeSens = 0x01
    typeMed = 0x10
    typeLrg = 0x11
    typeBrst = 0x20

    for obj in mypackets:
        if(obj.pcktType == typeSens):
            sensor.append(i)
        elif(obj.pcktType == typeMed):
            medium.append(i)
        elif(obj.pcktType == typeLrg):
            large.append(i)
        elif(obj.pcktType == typeBrst):
            burst.append(i)
        i=i+1
    return sensor, medium, large, burst

def ft_to_km(ft): #feet to kilometers
    return 1.609 * ft / 5280

def ascent_desent(h_func,t): #seperate out the ascent and descent parts of the flight
    h = []
    for i in t:
        h.append(ft_to_km(h_func(i/1000))) #function expects time in sec t in ms
    max_h = np.argmax(h) #get index of max height
    return t,h,max_h

def order_merge_3(time,data):#merges 3 arrays together in ascending time, keeps data with time
    n1 = len(time[0])
    n2 = len(time[1])
    n3 = len(time[2])
    ret_time = [None] * (n1 + n2 + n3)
    ret_data = [None] * (n1 + n2 + n3)
    i1 = 0
    i2 = 0
    i3 = 0
    k = 0
    while (i1 < n1 and i2 < n2 and i3 < n3):
        if (time[0][i1] > time[1][i2]):
            if (time[1][i2] > time[2][i3]):
                ret_time[k] = time[2][i3]
                ret_data[k] = data[2][i3]
                k = k + 1
                i3 = i3 + 1
            else:
                ret_time[k] = time[1][i2]
                ret_data[k] = data[1][i2]
                k = k + 1
                i2 = i2 + 1
        else:
            if (time[0][i1] > time[2][i3]):
                ret_time[k] = time[2][i3]
                ret_data[k] = data[2][i3]
                k = k + 1
                i3 = i3 + 1
            else:
                ret_time[k] = time[0][i1]
                ret_data[k] = data[0][i1]
                k = k + 1
                i1 = i1 + 1
    while (i1 < n1 and i2 < n2):
        if (time[0][i1] > time[1][i2]):
            ret_time[k] = time[1][i2]
            ret_data[k] = data[1][i2]
            k = k + 1
            i2 = i2 + 1
        else:
            ret_time[k] = time[0][i1]
            ret_data[k] = data[0][i1]
            k = k + 1
            i1 = i1 + 1
    while (i1 < n1 and i3 < n3):
        if (time[0][i1] > time[2][i3]):
            ret_time[k] = time[2][i3]
            ret_data[k] = data[2][i3]
            k = k + 1
            i3 = i3 + 1
        else:
            ret_time[k] = time[0][i1]
            ret_data[k] = data[0][i1]
            k = k + 1
            i1 = i1 + 1
    while (i3 < n3 and i2 < n2):
        if (time[2][i3] > time[1][i2]):
            ret_time[k] = time[1][i2]
            ret_data[k] = data[1][i2]
            k = k + 1
            i2 = i2 + 1
        else:
            ret_time[k] = time[2][i3]
            ret_data[k] = data[2][i3]
            k = k + 1
            i3 = i3 + 1
    while (i1 < n1):
        ret_time[k] = time[0][i1]
        ret_data[k] = data[0][i1]
        k = k + 1
        i1 = i1 + 1
    while (i2 < n2):
        ret_time[k] = time[1][i2]
        ret_data[k] = data[1][i2]
        k = k + 1
        i2 = i2 + 1
    while (i3 < n3):
        ret_time[k] = time[2][i3]
        ret_data[k] = data[2][i3]
        k = k + 1
        i3 = i3 + 1
    return ret_time,ret_data #returns sorted and merged list of data




def main():
    print('hello')
    return
