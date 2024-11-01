# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:52:56 2022

@author: calvin
"""

import numpy as np

#from struct import unpack

from scipy import optimize

from scipy import stats

#import pickle

import matplotlib.pyplot as plt

from matplotlib import gridspec

import UDIP_Lib_V20 as UDIP_Lib  #making updating UDIP_Lib easier

#from datetime import datetime

#1. Initialization 

filename = "UDIP0000.DAT"
Outpath = "Analysis/"

UDIP_packets = []
ind_sensor = []
ind_medium = []
ind_large = []
ind_burst = []
ind_constant = []


def Load_data():
    UDIP_read(filename)
    Find_Indexs()
    return

def UDIP_read(filename):
    global UDIP_packets 
    UDIP_packets = np.array(UDIP_Lib.readFile(filename))

#2. First Scan

def Plot_Precentiles(start=None,stop=None,percent=95):
    Packet_test()#load data only if needed
    #need to be updated!!!!!!!
    
    S_type = UDIP_Lib.typeMed
    
    time=[]
    adc0=[]
    adc1=[]
    adc2=[]
    
    if start == None and stop == None:
       for obj in UDIP_packets:
           if obj.pcktType == S_type:
               time.append(obj.tInitial)
               adc0.append(np.percentile(obj.sweep.adc0Curr, percent))
               adc1.append(np.percentile(obj.sweep.adc1Curr, percent))
               adc2.append(np.percentile(obj.sweep.adc2Curr, percent))
    
    ## this means if secStart and secStop are all set to 'None';
    ## then all objects will be added to the list and plotted
  
    elif start == None:
        for obj in UDIP_packets:
            if obj.pcktType == S_type and obj.tInitial <= 1000*stop:
                time.append(obj.tInitial)
                adc0.append(np.percentile(obj.sweep.adc0Curr, percent))
                adc1.append(np.percentile(obj.sweep.adc1Curr, percent))
                adc2.append(np.percentile(obj.sweep.adc2Curr, percent))
    
    ## this means if only secStart and is set to 'None';
    ## then the objects whose obj.tInitial attribute between 0 and secStop value
    ## will be added to the list and plotted
    
    elif stop == None:
        for obj in UDIP_packets:
            if obj.pcktType == S_type and obj.tInitial >=1000*start:
                time.append(obj.tInitial)
                adc0.append(np.percentile(obj.sweep.adc0Curr, percent))
                adc1.append(np.percentile(obj.sweep.adc1Curr, percent))
                adc2.append(np.percentile(obj.sweep.adc2Curr, percent))
    
    ## this means if only SecStop is set to 'None';
    ## then the objects whose obj.tInitial attribute are from secStart value to the end of the packet
    ## will be plotted 
    
    
    else:
        for obj in UDIP_packets:
            if obj.tInitial >= 1000*start and obj.pcktType == S_type:
                if obj.tInitial >= stop*1000:
                    break
                else:
                    time.append(obj.tInitial)
                    adc0.append(np.percentile(obj.sweep.adc0Curr, percent))
                    adc1.append(np.percentile(obj.sweep.adc1Curr, percent))
                    adc2.append(np.percentile(obj.sweep.adc2Curr, percent))
                
    ## this means if only SecStop is set to 'None';
    ## then the objects whose obj.tInitial attribute are from secStart value to the end of the packet
    ## will be plotted 
                
      
    fig=plt.figure(figsize=(10,10)) #left is base, right is height
    
    
    ax1 = fig.add_axes([0.20, 0.10, 0.6, 0.25])
    ax2 = fig.add_axes([0.20, 0.38, 0.6, 0.25])
    ax3 = fig.add_axes([0.20, 0.66, 0.6, 0.25])
    
    #plot the data 
    ax1.plot(time, adc0, 'firebrick', linewidth=1.6, label='ADC 0')
    ax2.plot(time, adc1, 'violet', linewidth=1.6, label='ADC 1')
    ax3.plot(time, adc2, 'olive', linewidth=1.6, label='ADC 2')
    
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper left")
    ax3.legend(loc="upper left")
    
    ax1.set_ylabel('Current (nA)', fontsize=12, color='blue')
    ax2.set_ylabel('Current (nA)', fontsize=12, color='blue')
    ax3.set_ylabel('Current (nA)', fontsize=12, color='blue')

    ax1.set_xlabel('Time (ms)', fontsize=12, color='blue')
    
    #plt.plot(time,adc0,'firebrick', linewidth=1.6, label='ADC 0')
    #plt.plot(time,adc1,'violet', linewidth=1.6, label='ADC 1')
    #plt.plot(time,adc2,'olive', linewidth=1.6, label='ADC 2')
    
    #plt.xlabel("Time")
    #plt.ylabel("Current")
    #plt.legend()
    
    if (S_type==UDIP_Lib.typeMed):
        typ = 'Medium'
    elif (S_type==UDIP_Lib.typeLrg):
        typ = 'Large'
    elif (S_type==UDIP_Lib.typeBrst):
        typ = 'Burst'
    fig.savefig(Outpath+'{}Sweep_from_{}sec_to_{}sec_{}_percentile.png'.format(typ,start,stop,percent))
    plt.close()


def Plot_Sensors(start,stop,accel=False,accelH=False,gyro=False,magno=False,temp=False,photo=False):
    Packet_test()#load data only if needed
    
    time_started = []
    target = []
    AccX = []
    AccY = []
    AccZ = []
    AccH = []
    #AccZ_ana = []
    GyroX = []
    GyroY = []
    GyroZ = []
    MagX = []
    MagY = []
    MagZ = []
    Temperature_p = []
    Temperature_s = []
    Temperature_d = []
    Pd_1 = []
    Pd_2 = []
    
    beggining = find_closest(start, ind_sensor)
    end = find_closest(stop, ind_sensor)
    
    for i in ind_sensor:
        if (i >= ind_sensor[beggining]) and (i <= ind_sensor[end]):
            time_started.append(UDIP_packets[i].tInitial)
            target.append(i)
        
    for i in target:
        AccX.append(UDIP_packets[i].accX)
        AccY.append(UDIP_packets[i].accY)
        AccZ.append(UDIP_packets[i].accZ)
        AccH.append(UDIP_packets[i].accH)
        #AccZ_ana.append(UDIP_packets[i].accAna)
        GyroX.append(UDIP_packets[i].gyroX)
        GyroY.append(UDIP_packets[i].gyroY)
        GyroZ.append(UDIP_packets[i].gyroZ)
        MagX.append(UDIP_packets[i].magX)
        MagY.append(UDIP_packets[i].magY)
        MagZ.append(UDIP_packets[i].magZ)
        Temperature_p.append(UDIP_packets[i].temperature_p)
        Temperature_s.append(UDIP_packets[i].temperature_s)
        Temperature_d.append(UDIP_packets[i].temperature_d)
        Pd_1.append(UDIP_packets[i].pd_1)
        Pd_2.append(UDIP_packets[i].pd_2)
    
    if (accel):
        plt.plot(np.divide(time_started,1000),AccX,label="X")
        plt.plot(np.divide(time_started,1000),AccY,label="Y")
        plt.plot(np.divide(time_started,1000),AccZ,label="Z")
        if (accelH):
            plt.plot(np.divide(time_started,1000),AccH,label="H")
        plt.title("Accelerometer")
        #plt.plot(np.divide(time_started,1000),AccZ_ana,label="Z Acceleration analog")
        plt.ylabel("Acceleration (m/s^2)")
        plt.xlabel("Time (sec)")
        plt.legend(loc="lower right")
        plt.savefig(Outpath+'acceleration_from_{}sec_to_{}sec.png'.format(start,stop))
        plt.close()
    
    if (gyro):
        plt.plot(np.divide(time_started,1000),GyroX,label="X")
        plt.plot(np.divide(time_started,1000),GyroY,label="Y")
        plt.plot(np.divide(time_started,1000),GyroZ,label="Z")
        plt.title('Gyroscope')
        plt.ylabel("Spin rate (degrees/sec)")
        plt.xlabel("Time (sec)")
        plt.legend(loc="upper left")
        plt.savefig(Outpath+'gyro_from_{}sec_to_{}sec.png'.format(start,stop))
        plt.close()
    
    if (magno):
        plt.plot(np.divide(time_started,1000),MagX,label="X")
        plt.plot(np.divide(time_started,1000),MagY,label="Y")
        plt.plot(np.divide(time_started,1000),MagZ,label="Z")
        plt.title('Magnetometer')
        plt.ylabel("Magnetic Field (Guass)")
        plt.xlabel("Time (sec)")
        plt.legend(loc="upper left")
        plt.savefig(Outpath+'mag_from_{}sec_to_{}sec.png'.format(start,stop))
        plt.close()
    
    if (temp):
        plt.plot(np.divide(time_started,1000),Temperature_p,label="Pocket")
        plt.plot(np.divide(time_started,1000),Temperature_s,label="Skin")
        plt.plot(np.divide(time_started,1000),Temperature_d,label="Digital")
        plt.title('Temperature')
        plt.ylabel("Temperature (C)")
        plt.xlabel("Time (sec)")
        plt.legend(loc="upper left")
        plt.savefig(Outpath+'temperature_from_{}sec_to_{}sec.png'.format(start,stop))
        plt.close()
    
    if (photo):
        plt.plot(np.divide(time_started,1000),Pd_1,label="Pd 1")
        plt.plot(np.divide(time_started,1000),Pd_2,label="Pd 2")
        plt.title('Photodiode')
        plt.ylabel("???")
        plt.xlabel("Time (sec)")
        plt.legend(loc="upper left")
        plt.savefig(Outpath+'photo_diodes_from_{}sec_to_{}sec.png'.format(start,stop))
        plt.close()

#3. Packet Exploration

def Plot_sweep(time=None,index=None,fit=False,hysteresis=False):
    Packet_test()#load data only if needed
    if (time != None) or (index != None):
        if (index == None):
            ind = find_closest(time, ind_medium)
        else:
            ind = index
        ind = ind_medium[ind] #get index of medium packet
        packet = UDIP_packets[ind] #get medium packet
        
        #adc0
        x0 = packet.sweep.sweepVoltage
        y0 = packet.sweep.adc0Curr
        
        #adc1
        x1 = packet.sweep.sweepVoltage
        y1 = packet.sweep.adc1Curr
        
        #adc2
        x2 = packet.sweep.sweepVoltage
        y2 = packet.sweep.adc2Curr
        
        #photo diodes
        #i_pd_1 = packet.sweep.i_pd_1
        #i_pd_2 = packet.sweep.i_pd_2
        #f_pd_1 = packet.sweep.f_pd_1
        #f_pd_2 = packet.sweep.f_pd_2
        
        fig=plt.figure(figsize=(10,10)) 
        
        ax1 = fig.add_axes([0.3, 0.09, 0.6, 0.25])
        ax2 = fig.add_axes([0.3, 0.37, 0.6, 0.25])
        ax3 = fig.add_axes([0.3, 0.65, 0.6, 0.25])
        
        ax1.plot(x0, y0, 'firebrick', linewidth=1.6)#, label='Low Gain')
        if (fit):
            t, m ,popt,pcov = fit_sweep(x0, y0, hysteresis)
            ax1.plot(t,m,label='Model',linewidth=3)
        ax2.plot(x1,y1, 'violet', linewidth=1.6)#, label='Mid Gain')
        ax3.plot(x2,y2, 'olive', linewidth=1.6)#, label='High Gain')
        
        #ax1.legend(loc="upper left")
        #ax2.legend(loc="upper left")
        #ax3.legend(loc="upper left")
        
        ax1.set_ylabel(r'Probe Current, $I_{\rm probe}$ (nA)' '\n' r'Low Gain', fontsize=12, color='black')
        ax2.set_ylabel(r'Probe Current, $I_{\rm probe}$ (nA)' '\n' r'Mid Gain', fontsize=12, color='black')
        ax3.set_ylabel(r'Probe Current, $I_{\rm probe}$ (nA)' '\n' r'High Gain', fontsize=12, color='black')
        ax1.set_xlabel(r'Sweep Voltage, $V_{\rm sweep} (V)$', fontsize=12, color='black')
        
        #dashed line at x=0
        ax1.axvline(x= 0, color ="lightgray", linewidth=1.7, linestyle =":")
        ax2.axvline(x= 0, color ="lightgray", linewidth=1.7, linestyle =":")
        ax3.axvline(x= 0, color ="lightgray", linewidth=1.7, linestyle =":")
        
        #dashed line at y=0
        ax1.axhline(y= 0, color ="lightgray", linewidth=1.7, linestyle =":")
        ax2.axhline(y= 0, color ="lightgray", linewidth=1.7, linestyle =":")
        ax3.axhline(y= 0, color ="lightgray", linewidth=1.7, linestyle =":")
        plt.show()
        plt.close()
        
def Plot_Many_Sweeps():
    Packet_test()#load data only if needed
    print("The number of medium sweeps is " + str(len(ind_medium)))
    while(True):
        cont = input("enter q to end the loop")
        if (cont == "q"):
            break
        ind = int(input("enter the integer index of the medium sweep you want to see "))
        f = input("enter t if you want the data to be fit")
        if (f == "t"):
            fit = True
        else:
            fit = False
        h = input("enter t if you want hysteresis to be removed")
        if (h == "t"):
            hyst = True
        else:
            hyst = False
        Plot_sweep(time=None,index=ind,fit=fit,hysteresis=hyst)

def Plot_constant(time=None,index=None):
    Packet_test()#load data only if needed
    if (time != None) or (index != None):
        if (index == None):
            ind = find_closest(time, ind_constant)
        else:
            ind = index
        ind = ind_constant[ind] #get index of constant packet
        packet = UDIP_packets[ind] #get constant packet
        
        #photodiode
        pd1_val = packet.sweep.adcPD1Val
        pd2_val = packet.sweep.adcPD2Val
        
        #voltages
        pd1_v = packet.sweep.adcPD1Volt
        pd2_v = packet.sweep.adcPD2Volt

        adc0_c = packet.sweep.adc0Curr
        
        #print(np.average(packet.sweep.ProbeVolt))

        # fig = plt.figure()
        # gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1]) 
        # ax0 = plt.subplot(gs[0])
        # line0, = ax0.plot(pd1_v,label="pd1")
        # line0, = ax0.plot(pd2_v,label="pd2")

        # ax1 = plt.subplot(gs[1], sharex = ax0)
        # line1, = ax1.plot(adc0_c,label="adc0_c")

        fig, axs = plt.subplots(2)
        fig.suptitle(f'Photodiode vs Current\nConstant Voltage: {str(np.round(np.average(packet.sweep.ProbeVolt),3))} V\nTime (ms): {packet.tInitial}')
        axs[0].plot(pd1_v,label="pd1")
        axs[0].plot(pd2_v,label="pd2")
        axs[0].legend()
        axs[0].set_ylabel('Photodiode (V)')
        axs[1].plot(adc0_c,label="adc0_c")
        axs[1].legend()
        axs[1].set_ylabel('Current (nA)')

        
        # plt.plot(pd1_val,pd1_v,label="pd1")
        # plt.plot(pd2_val,pd2_v,label="pd2")
        # plt.plot(pd1_v,label="pd1")
        # plt.plot(pd2_v,label="pd2")
        #plt.title("Constant Voltage " + str(np.round(np.average(packet.sweep.ProbeVolt),3)) + "V")
        #plt.suptitle(f"Time (ms): {packet.tInitial}")
        #ax0.ylabel("Photodiodes (V)")
        #ax1.ylabel("ADC0 Current (nA)")
        plt.legend()
        plt.subplots_adjust(hspace=.0)
        plt.savefig(f'{Outpath}const_sweep_from_{packet.tInitial}sec.png')
        #plt.show()
        plt.close()

#4. helper functions

def Remove_hysteresis(x,y):
    print("to be completed")
    return x,y

def Find_Indexs():
    global ind_sensor
    global ind_medium
    global ind_large
    global ind_burst
    global ind_constant
    i = 0
    #pcktType
    typeSens = 0x01
    typeMed = 0x10
    typeLrg = 0x11
    typeBrst = 0x20
    typeCst = 0x30
    
    for obj in UDIP_packets:
        if(obj.pcktType == typeSens):
            ind_sensor.append(i)
        elif(obj.pcktType == typeMed):
            ind_medium.append(i)
        elif(obj.pcktType == typeLrg):
            ind_large.append(i)
        elif(obj.pcktType == typeBrst):
            ind_burst.append(i)
        elif(obj.pcktType == typeCst):
            ind_constant.append(i)
        i=i+1

def Packet_test():
    #test to see if UDIP_packets is loaded
    if (UDIP_packets == []):
        Load_data()

def model(x,xa,b,m1,n,t,V0):
    #electron retardation
    ret = np.zeros(len(x))
    ret[x <= xa] = seg1(x[x <= xa],m1) - seg1(xa, m1) + b
    ret[x > xa] = seg2(x[x > xa],n,t,V0) - seg2(xa,n,t,V0) + seg1(xa,m1) - seg1(xa, m1) + b
    return ret
    
def seg1(x,m):#linear--full model square root
    return m * x

def seg2(x,n,t,V0):# square root
    q_e = 1.602 * 10**-19 #C                charge of an electron
    K_b = 1.381 * 10**-23 #m^2*kg/(s^2*K)   boltzman constant
    m_e = 9.109 * 10**-31 #kg               mass of an electron
    R = (3./16.) * 0.0254 #radius of probe
    L = (3.25) * 0.0254 #length of probe
    A = 2. * np.pi * R * L + np.pi * (R ** 2) #area of probe cylinder with out a bottom
    
    k = q_e / (K_b * t)
    I0 =n * q_e * np.sqrt(K_b * t / (2. * np.pi * m_e)) * A / (10**-9)
    return I0 * np.sqrt(1. + k*(x + V0))

def get_fit(x,y):
    g = [0.6,-14,80, 5*(10**10),1000,-0.5]    #intial guess
    b = ((-3,-np.inf,-np.inf,0,0,-3),(3,np.inf,np.inf,np.inf,10000,3)) #bounds
    popt, pcov = optimize.curve_fit(model,x,y,g,bounds=b)
    #print(popt)
    max_1 = max(x)
    min_1 = min(x)
    t = np.linspace(min_1,max_1,num=60)
    return t, model(t,*popt),popt,pcov #popt[0:xa,1:b,2:m1,3:n,4:t,5:V0]

def noise_cleaning(x_raw,y_raw):
    #data filtering goes here
    return x_raw,y_raw

def fit_sweep(x_raw,y_raw,hysteresis):
    x,y = noise_cleaning(x_raw, y_raw)
    if (hysteresis):
        x,y = Remove_hysteresis(x, y)
    return get_fit(x, y)

def find_closest(time,array):
    #returns index in array that is closest to the given time
    n = len(array)
    if (time < UDIP_packets[array[0]].tInitial):
        return 0
    elif (time > UDIP_packets[array[-1]].tInitial):
        return n-1
    else:
        b = UDIP_packets[array[0]].tInitial
        a = UDIP_packets[array[0]].tInitial
        i = 1
        while (a < time):
            b = a
            a = UDIP_packets[array[i + 1]].tInitial
            i = i + 1
        if (time - b < a - time):
            return i - 1
        else:
            return i

#Plot_sweep(index=500)
#Plot_Many_Sweeps()

Plot_Precentiles()

Plot_Sensors(0, 300000,accel=True,accelH=True,gyro=True,magno=True,temp=True,photo=True)
Plot_Sensors(0, 50000,accel=True,accelH=True,gyro=True,magno=True)
Plot_Sensors(45000, 46000,gyro=True,magno=True,photo=True)
Plot_Sensors(0, 600000,temp=True)
Plot_Sensors(0, 1000000,temp=True)
Plot_Sensors(0, 8000000,temp=True)
Plot_Sensors(600000, 800000,temp=True)
Plot_Sensors(0, 4000,accel=True)
Plot_Sensors(700000, 1000000,accel=True)

Plot_constant(index=35)
Plot_constant(index=44)
Plot_constant(index=53)
Plot_constant(index=62)
Plot_constant(index=71)
Plot_constant(index=80)
Plot_constant(index=89)


