# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 20:35:53 2020
Class definitions for Reading Packets from University of Delaware Ionosphere
Probe RockSat-C 2022 Mission
Last updated 22/05/23
@author: Jarod Dagney
"""
from numpy import tile

from struct import unpack

from scipy import interpolate

import math

#from matplotlib import pyplot as pt

import os

import matplotlib.pyplot as plt
#from scipy.io import savemat
#import numpy as np

VERBOSE = False
lenHedr = 15

lenSens = 30
lenBrst = 133 * 8 + 12
lenMed = 265 * 8 + 12
lenLrg = 253 * 8 + 12
lenCst = 100 * 10 + 6

typeSens = 0x01
typeMed  = 0x10
typeLrg  = 0x11
typeBrst = 0x20
typeCst  = 0x30
#------------------------------------------------------#
#Parent Class of all Packets
#Attributes:
#Global:
#   totCnt- Total number of created packets
#Local:
#   tInitial- Time of packet creation in ms (float)
#   tFinal- Final time of packet creation in ms (float)
#   pcktType- Type of Packet (int)
#       Sensor = 0x01        Medium Sweep = 0x10
#       Large Sweep = 0x11   Burst Sweep = 0x20
#   payloadLen- Length of packet payload (int)
#------------------------------------------------------#
class Packet:
    """
    Parent Class of all Packets. Subclasses: Sensor_Packet, Medium_Sweep,
    Large_Sweep, Burst_Sweep.
    
    Attributes:
        Global:
            totCnt- Total number of created packets
        Local:
            count- This object's count
            tInitial- Time of packet creation in ms (float)
            tFinal- Final time of packet creation in ms (float)
            pcktType- Type of Packet (int)
                Sensor = 0x01        Medium Sweep = 0x10
                Large Sweep = 0x11   Burst Sweep = 0x20
            payloadLen- Lenght of packet payload (int)
    """
    totCnt = 0
    def __init__(self, count, tInitial, tFinal, myType, payloadLen):
        '''
        Initialization function for class: Packet.
        
        Parameters:
            count: 2 byte unsigned int 
            tInitial: 4 byte unsigned int
            tFinal: 4 byte unsigned int
            myType: 1 byte "char" for redundancy on packet types
            payloadLen: 2 byte int, total number of bytes in the payload
        '''
        Packet.totCnt += 1
        self.count = count
        self.tInitial = tInitial
        self.tFinal = tFinal
        self.pcktType = myType
        self.payloadLen = payloadLen
        
    
    def __del__(self):
        type(self).totCnt -= 1
        
    
    

class Sensor_Packet(Packet):
    """
    Subclass of Packet, holds the sensor packet data.
    
    Attributes:
        All Packet attributes
        pyld: List of bytes to convert
        accelScale: Maximum range of accelerometer 2, 4, 8, or 16 g
        accX, accY, accZ: Acceleration in g (9.81 m/s^2) in three axes
        gyrX, gyrY, gyrZ: Spin rate in degrees per second in three axes
        magX, magY, magZ: Magnetic field in Gauss in three axes
        accAna: Accerleration in g (9.81 m/s^2) in z axis (along path of rocket) from analog sensor
        temperature: Internal canister temperature
        temperatureAna: Internal canister temperature from analog sensor
    """
    SENSITIVITY_ACCELEROMETER_2 = 0.000061
    SENSITIVITY_ACCELEROMETER_4 = 0.000122
    SENSITIVITY_ACCELEROMETER_8 = 0.000244
    SENSITIVITY_ACCELEROMETER_16= 0.000732
    
    SENSITIVITY_GYROSCOPE_2000  = 0.07
    
    SENSITIVITY_MAGNETOMETER_4  = 0.00014
    
    temperatureConversion = 0.0625

    def __init__(self, count, tInitial, tFinal, pcktType, payloadLen, pyld):
        """
        Initialization method for class: Sensor_Packet.
        
        Parameters:
            count: 2 byte unsigned int 
            tInitial: 4 byte unsigned int
            tFinal: 4 byte unsigned int
            myType: 1 byte "char" for redundancy on packet types
            payloadLen: 2 byte int, total number of bytes in the payload
            pyld: List of bytes
        """
        super().__init__(count, tInitial, tFinal, pcktType, payloadLen)
        self.pyld = pyld
        self.readPyld()
        
           
    
    def readPyld(self):
        """Extracts bytes from pyld and populates sensor attributes"""
        #Unpack Acceleration data from payload
        #self.accScale = unpack('<B', self.pyld[lenSens - 1:lenSens])[0]
        #self.acc = unpack('<hhhH',self.pyld[0:8])#three signed two byte integers and 1 unsigned two byte integer
        self.acc_m = unpack('<hhh',self.pyld[0:6])
        self.acc_h = unpack("<h", self.pyld[6:8])[0]
        self.calibrateDigAcc()

        #Acceleration data from analog sensor
        
        
        #Unpack Gyroscope data from payload
        self.gyr = unpack('<hhh', self.pyld[8:14])
        self.calibrateDigGyr()
        
        #Unpack Magnetic Field data from payload
        self.mag = unpack('<hhh', self.pyld[14:20])
        self.calibrateDigMag()
        
        #Unpack Temperature data from payload
        tmp_d, tmp_p, tmp_s = unpack('<HHH', self.pyld[20:26])
        self.temperature_d = tmp_d * self.temperatureConversion
        b = 3455
        if(tmp_p == 0):
            print('Here')
            tmp_p = 0.1
        if(tmp_s == 0):
            tmp_s = 0.1
        tmp_p_res = ((3.3 * 27000) / ((tmp_p*3.3)/(4096*16))) - 27000
        tmp_s_res = ((3.3 * 27000) / ((tmp_s*3.3)/(4096*16))) - 27000
        self.temperature_p = b / math.log(tmp_p_res / (10000 * math.exp((-1 * b) / 298 ))) - 273
        self.temperature_s = b / math.log(tmp_s_res / (10000 * math.exp((-1 * b) / 298 ))) - 273
        
        #Unpack photodiode data
        self.pd_1, self.pd_2 = unpack('<HH', self.pyld[26:30])
        
    def calibrateDigAcc(self):
        self.accScale = 16 
        """Calibrates the acceleration and stores it in respective class attributes in g"""
        self.accH = self.acc_h / 33.4
        # if(self.accScale == 2):
        #     self.accX = self.acc_m[0] * self.SENSITIVITY_ACCELEROMETER_2
        #     self.accY = self.acc_m[1] * self.SENSITIVITY_ACCELEROMETER_2
        #     self.accZ = self.acc_m[2] * self.SENSITIVITY_ACCELEROMETER_2
        # elif(self.accScale == 4):
        #     self.accX = self.acc_m[0] * self.SENSITIVITY_ACCELEROMETER_4
        #     self.accY = self.acc_m[1] * self.SENSITIVITY_ACCELEROMETER_4
        #     self.accZ = self.acc_m[2] * self.SENSITIVITY_ACCELEROMETER_4
        # elif(self.accScale == 8):
        #     self.accX = self.acc_m[0] * self.SENSITIVITY_ACCELEROMETER_8
        #     self.accY = self.acc_m[1] * self.SENSITIVITY_ACCELEROMETER_8
        #     self.accZ = self.acc_m[2] * self.SENSITIVITY_ACCELEROMETER_8
        # elif(self.accScale == 16):
        #     self.accX = self.acc_m[0] * self.SENSITIVITY_ACCELEROMETER_16
        #     self.accY = self.acc_m[1] * self.SENSITIVITY_ACCELEROMETER_16
        #     self.accZ = self.acc_m[2] * self.SENSITIVITY_ACCELEROMETER_16
        # else:
        #     #Not valid accel scale
        #     self.accX = 0
        #     self.accY = 0
        #     self.accZ = 0
        self.accX = self.acc_m[0] / 95.43
        self.accY = self.acc_m[1] / 95.43
        self.accZ = self.acc_m[2] / 95.43
            
    
    def calibrateDigGyr(self):
        """Calibrates the spin rate and stores it in respective class attributes in degrees per second"""
        self.gyroX = self.gyr[0] / 936.25
        self.gyroY = self.gyr[1] / 936.25
        self.gyroZ = self.gyr[2] / 936.25
        
    
    def calibrateDigMag(self):
        """Calibrates the magnetic field data and stores it in respective class attributes in Gauss"""
        self.magX = self.mag[0] / 409.6
        self.magY = self.mag[1] / 409.6
        self.magZ = self.mag[2] / 409.6
        
    

class Large_Sweep(Packet):
    """
    Subclass of Packet, creates an object of type Sweep.
    
    Attributes:
        All Packet Attributes
        sweep: Object of type Sweep that holds Voltage and Current values
    """
    nSteps = 253
    
    def __init__(self, count, tInitial, tFinal, pcktType, payloadLen, pyld):
        """
        Initialization method for class: Large_Sweep
        
        Parameters:
            count: 2 byte unsigned int 
            tInitial: 4 byte unsigned int
            tFinal: 4 byte unsigned int
            myType: 1 byte "char" for redundancy on packet types
            payloadLen: 2 byte int, total number of bytes in the payload
            pyld: List of bytes 
        """
        super().__init__(count, tInitial,tFinal,pcktType,payloadLen)
        self.sweep = Sweep(pyld,self.nSteps)
        return

class Medium_Sweep(Packet):
    """
    Subclass of Packet, creates an object of type Sweep.
    
    Attributes:
        All Packet Attributes
        sweep: Object of type Sweep that holds Voltage and Current values
    """
    
    nSteps = 265
    
    def __init__(self, count, tInitial, tFinal, pcktType, payloadLen, pyld):
        super().__init__(count, tInitial, tFinal, pcktType, payloadLen)
        self.sweep = Sweep(pyld,self.nSteps)
        return

class Burst_Sweep(Packet):
    """
    Subclass of Packet, creates an array of objects of type Sweep.
    
    Attributes:
        All Packet Attributes
        pyld
        sweepArray: List of objects of type Sweep that holds Voltage and Current values
    """
    nSteps = 133
    
    def __init__(self, count, tInitial, tFinal, pcktType, payloadLen, pyld):
        super().__init__(count, tInitial, tFinal, pcktType, payloadLen)
        vRef = pyld[0:4]
        self.sweepArray = []
        for i in range(10):
            self.sweepArray.append(Sweep(vRef + pyld[ 4 + (i*133*2*4) : 4+(i*133*2*4)+(133*2*4) ], self.nSteps))
        
        return

class Constant_Sweep(Packet):
    """
    Subclass of Packet, creates an object of type Sweep.
    
    Attributes:
        All Packet Attributes
        sweep: Object of type Sweep that holds Voltage and Current values
    """
    
    nSteps = 100
    
    def __init__(self, count, tInitial, tFinal, pcktType, payloadLen, pyld):
        super().__init__(count, tInitial, tFinal, pcktType, payloadLen)
        self.sweep = CstSweep(pyld,self.nSteps)
        return
       

class Sweep():
    """
    Holds Voltage and Current values from a voltage sweep.
    
    Attributes:
        R* are resistor values on board analog board
        pyld: List of bytes to convert
        nSteps: Number of steps in the voltage sweep
        vPos, vNeg: Voltage supply to analog board used as reference
        sweepVoltage: List of voltages applied to the probe
        adc0Curr: List of currents induced by the probe
        adc1Curr: List of currents induced by the probe with one gain stage
        adc2Curr: List of currents induced by the probe with two gain stages
        
    """
    nStepsMed = 265
    nStepsLrg = 253
    nStepsBrst = 133 * 10
    
    nAvg = 16.0
    
    R1 = 10000.0
    R2 = 10000.0
    R3 = 30000.0
    R4 = 10000.0
    
    RG1 = 10000.0
    RG2 = 191000.0
    
    RS1 = 33200
    RS2 = 4990
    RS3 = 100000
    
    ###For Sweep Voltage Divider, NOT in use as of 3/18/21
    #RS1 = 30100.0
    #RS2 = 49900.0
    
    #vPos Voltage Divider Resistors
    RP1 = 49900.0
    RP2 = 10000.0
    
    #vNeg Voltage Divider Resistors
    RN1 = 49900.0
    RN2 = 69800.0
    
    RF = 1000000.0
    
    def __init__(self, pyld, nSteps):
        """
        Initialization method of class Sweep.
        
        Parameters:
            pyld: List of bytes
            nSteps: Number of steps in the voltage sweep. Sets up the length of
            attributes
        """
        self.nSteps = nSteps
        self.pyld = pyld
        
        self.vPosVal = unpack('<H', pyld[0:2])[0]
        self.vNegVal = unpack('<H', pyld[2:4])[0]

        self.i_pd_1, self.i_pd_2 = unpack('<HH', pyld[4:8])
        self.f_pd_1, self.f_pd_2 = unpack('<HH', pyld[8:12])
        
        self.vPos = (self.vPosVal / 16.0) * (3.3 / 4096.0) * (self.RP1 + self.RP2) / self.RP2
        self.vNeg = (self.RN1 + self.RN2) / self.RN1 * ( ( (self.vNegVal / 16.0) * (3.3 / 4096.0) ) - (self.vPos ) * (self.RN2 / (self.RN1 + self.RN2)))
        
        self.adcDACVal =  tile(0. , nSteps)
        self.adcDACVolt = tile(0., nSteps)
        self.sweepVoltage = tile(0., nSteps)
        
        self.adc0Val =  tile(0. , nSteps)
        self.adc0Volt = tile(0., nSteps)
        self.adc0Curr = tile(0., nSteps)
        
        self.adc1Val =  tile(0. , nSteps)
        self.adc1Volt = tile(0., nSteps)
        self.adc1Curr = tile(0., nSteps)
        
        self.adc2Val =  tile(0. , nSteps)
        self.adc2Volt = tile(0., nSteps)
        self.adc2Curr = tile(0., nSteps)

        self.ProbeVolt = tile(0. , self.nSteps)
        
        self.adc0CurrMax = 0
        self.adc1CurrMax = 0
        self.adc2CurrMax = 0
        self.adc0CurrMin = 0
        self.adc1CurrMin = 0
        self.adc2CurrMin = 0
        
        self.fillVals()
        self.fillVolt()
        self.fillCurr()
        self.sweepTransfer()
        return
            
    
    def setDACcmd(dacArr):
        x = 1
        return
    
    def fillVals(self):
        """Unpacks the DAC voltage and three ADC values"""
        for i in range(self.nSteps):
            arr = unpack('<HHHH', self.pyld[(4*2*i) + 12: (4*2*i + 8) + 12])
            self.adcDACVal[i] = arr[0] / self.nAvg
            #print(arr[0])
            self.adc0Val[i] = arr[1] / self.nAvg
            #print(arr[1])
            self.adc1Val[i] = arr[2] / self.nAvg
            self.adc2Val[i] = arr[3] / self.nAvg
            
        return

    
    def fillVolt(self):
        """Converts a value from ADC to the voltage"""
        conv = 3.3 / 4095.
        for i in range(self.nSteps):
            self.adcDACVolt[i] = self.adcDACVal[i] * conv
            self.adc0Volt[i] = self.adc0Val[i] * conv
            self.adc1Volt[i] = self.adc1Val[i] * conv
            self.adc2Volt[i] = self.adc2Val[i] * conv
            self.ProbeVolt[i] = (self.adcDACVal[i] * conv - 1.5) * 24
        return
         
    def fillCurr(self):
        """"Converts the voltage to ADC to a current from the probe"""
        for i in range(self.nSteps):
           #self.adc0Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc0Volt[i] - (self.R1/(self.R1+self.R2)) * 11.97 ) * ((self.RG1/(self.RG1+self.RG2))**0) * 1E9
           #self.adc1Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc1Volt[i] - (self.R1/(self.R1+self.R2)) * 11.97 ) * ((self.RG1/(self.RG1+self.RG2))**1) * 1E9
           #self.adc2Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc2Volt[i] - (self.R1/(self.R1+self.R2)) * 11.97 ) * ((self.RG1/(self.RG1+self.RG2))**2) * 1E9
            self.adc0Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc0Volt[i] - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**0) * 1E9
            self.adc1Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc1Volt[i] - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**1) * 1E9
            self.adc2Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc2Volt[i] - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**2) * 1E9
        self.adc0CurrMax = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 3.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**0) * 1E9
        self.adc1CurrMax = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 3.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**1) * 1E9
        self.adc2CurrMax = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 3.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**2) * 1E9
        self.adc0CurrMin = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 0.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**0) * 1E9
        self.adc1CurrMin = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 0.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**1) * 1E9
        self.adc2CurrMin = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 0.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**2) * 1E9
        
        return
           
        
    def sweepTransfer(self):
        """Converts voltage to ADC from DAC pin to voltage applied to probe"""
        #self.sweepVoltage = (self.RS3/self.RS1 + self.RS3/self.RS2 + 1) * self.adcSweepVolt - (self.RS3/self.RS1) * 11.97
        #self.sweepVoltage = ( ( self.RS1 + self.RS2 ) / self.RS2 * self.adcSweepVolt) - ( self.RS1 / self.RS2 ) * self.vNeg
        self.sweepVoltage = (self.RS3/self.RS1 + self.RS3/self.RS2 + 1) * self.adcDACVolt - (self.RS3/self.RS1) * self.vPos
        return
    
class CstSweep(Sweep):

    def fillVals(self):
        """Unpacks the DAC voltage and three ADC values"""
        #self.adcDACVal = unpack('<H', self.pyld[4:6])*self.nSteps
        self.adcDACVal = [(x) / self.nAvg for x in unpack('<H', self.pyld[4:6])*self.nSteps]
        self.adcPD1Val = tile(0. , self.nSteps)
        self.adcPD2Val = tile(0. , self.nSteps)
        self.adcPD1Volt = tile(0. , self.nSteps)
        self.adcPD2Volt = tile(0. , self.nSteps)
        self.ProbeVolt = tile(0. , self.nSteps)

        for i in range(self.nSteps):
            arr = unpack('<HHHHH', self.pyld[(5*2*i) + 6: (5*2*i + 10) + 6])
            #self.adcDACVal[i] = arr[0] / self.nAvg
            #print(arr[0])
            self.adcPD1Val[i] = arr[0] / self.nAvg
            self.adcPD2Val[i] = arr[1] / self.nAvg
            self.adc0Val[i] = arr[2] / self.nAvg
            #print(arr[1])
            self.adc1Val[i] = arr[3] / self.nAvg
            self.adc2Val[i] = arr[4] / self.nAvg
            
        return

    def fillVolt(self):
        """Converts a value from ADC to the voltage"""
        conv = 3.3 / 4095.
        for i in range(self.nSteps):
            self.adcDACVolt[i] = self.adcDACVal[i] * conv
            self.adc0Volt[i] = self.adc0Val[i] * conv
            self.adc1Volt[i] = self.adc1Val[i] * conv
            self.adc2Volt[i] = self.adc2Val[i] * conv
            self.adcPD1Volt[i] = self.adcPD1Val[i] * conv
            self.adcPD2Volt[i] = self.adcPD2Val[i] * conv
            self.ProbeVolt[i] = (self.adcDACVal[i] * conv - 1.5) * 24
        return

def readFile(fileName):
    """
    Main function of Conversion Software. Reads a file and creates a list of packet objects
    
    Parameters:
        fileName (str): Name of file
        
    Returns:
        myPackets (list): List of packet objects
    
    """
    myPackets = []
    myFile = open(fileName,"rb")
    raw = myFile.read()
    loc = 0
    i = 0
    while(loc < len(raw)):
        sync     = unpack('<BB',raw[loc + 0:loc + 2])
        count    = unpack('<H', raw[loc + 2:loc + 4])[0]
        tInitial = unpack('<I', raw[loc + 4:loc + 8])[0]
        tFinal   = unpack('<I', raw[loc + 8:loc + 12])[0]
        pcktType = unpack('<B', raw[loc + 12:loc + 13])[0]
        pyldLen  = unpack('<H', raw[loc + 13:loc + 15])[0]
        if(verifyHeader(sync, pcktType, pyldLen)):
            if(pcktType == typeSens):
                packet = Sensor_Packet(count, tInitial, tFinal, pcktType, pyldLen, raw[(loc+lenHedr) : (loc+lenHedr+pyldLen)] )
                myPackets.append(packet)
            elif(pcktType == typeMed):
                packet = Medium_Sweep(count, tInitial, tFinal, pcktType, pyldLen, raw[(loc+lenHedr) : (loc+lenHedr+pyldLen)])
                myPackets.append(packet)
            elif(pcktType == typeLrg):
                packet = Large_Sweep(count, tInitial, tFinal, pcktType, pyldLen, raw[(loc+lenHedr) : (loc+lenHedr+pyldLen)])
                myPackets.append(packet)
            elif(pcktType == typeBrst):
                packet = Burst_Sweep(count, tInitial, tFinal, pcktType, pyldLen, raw[(loc+lenHedr) : (loc+lenHedr+pyldLen)])
                myPackets.append(packet)
            elif(pcktType == typeCst):
                packet = Constant_Sweep(count, tInitial, tFinal, pcktType, pyldLen, raw[(loc+lenHedr) : (loc+lenHedr+pyldLen)])
                myPackets.append(packet)
            else:
                loc+=1
                print('Missing link')
            loc += lenHedr+pyldLen
            if(VERBOSE):
                print("Packet #: ", i)
            i+=1
        else:
            loc+=1
        
    return myPackets
                
    
def verifyHeader(sync, pcktType, pyldLen):
    """
    Helper function for readFile
    
    Parameters:
        sync (list): two byte array
        pcktType (int): type of packet in numerical value
        pyldLen (int): length of payload in bytes
        
    Returns:
        Boolean for if the header is valid
    """
    if(sync[0] != 0x55 or sync[1] != 0x44):
        return False
    if(pcktType == typeSens):
        if(pyldLen != lenSens):
            return False
        else:
            return True
    if(pcktType == typeMed):
        if(pyldLen != lenMed):
            return False
        else:
            return True
    if(pcktType == typeLrg):
        if(pyldLen != lenLrg):
            return False
        else:
            return True
    if(pcktType == typeBrst):
        if(pyldLen != lenBrst):
            return False
        else:
            return True
    if(pcktType == typeCst):
        if(pyldLen != lenCst):
            return False
        else:
            return True


def readTele(fileName):
    """
    Extracts the data from telemetry file
    
    Parameters
    ----------
    fileName : String
        Name of the file, should be telemetry file DQCA_quicklook...

    Returns
    -------
    flightTime : list
        list of the times from telemetry
    altitude : list
    
    horRange : list
    
    velocity : list

    """
    teleFile = open(fileName, "r")
    #get rid of the first 5 lines which hold no data
    for line in range(5):
        teleFile.readline()
    
    fileArr = teleFile.readlines()
    
    arrLen = len(fileArr)
    
    flightTime = []
    altitude = []
    horRange = []
    velocity = []
    
    for i in range(arrLen):
        data = fileArr[i].split()
        flightTime.append(float(data[0]))
        altitude.append(float(data[2]))
        horRange.append(float(data[3]))
        velocity.append(float(data[4]))
    
    return flightTime, altitude, horRange, velocity

def fitFlight(fileName):
    """
    Creates the fit function for the altitude as a function of time (seconds)

    Parameters
    ----------
    fileName : String
        DESCRIPTION.

    Returns
    -------
    Interpolatant function which takes a time (in seconds) and returns the altitude

    """
    flightTime, alt, rng, vel = readTele(fileName)
    
    return interpolate.interp1d(flightTime,alt), interpolate.interp1d(flightTime, vel)
    
#medSweep = readFile("medswp2.DAT")
#
#
#x = []
#for i in range(medSweep[0].sweep.nSteps-1):
#    x.append(medSweep[0].sweep.sweepVoltage[i]-medSweep[0].sweep.sweepVoltage[i+1])
#    
#plt.plot(range(264),x)
    


#for i in range(30):
    #medSweep[i].sweep.sweepVoltage


###mediumSweepFile = "MAKEMSWP.DAT"
###largeSweepFile = "MAKELSWP.DAT"
###burstSweepFile = "MAKEBSWP.DAT"
###sensorFile = "MAKESENS.DAT"

###mediumSweeps = readFile(mediumSweepFile)
###largeSweeps = readFile(largeSweepFile)
###burstSweeps = readFile(burstSweepFile)
###sensorPackets = readFile(sensorFile)


### mediumSweeps will be a list of 30 or 31 Medium_Sweep objects
### largeSweeps will be a list of 30 or 31 Large_Sweep objects
### and similarly for burstSweeps and sensorPackets
### To access one object in the list you can do mediumSweeps[i]
### or you can do a for each loop. Example: will print out the list of sweepVoltage of each Medium_
### Sweep object

#for obj in mediumSweeps:
#    print(obj.sweep.sweepVoltage)

#print(mediumSweeps[0].sweep.sweepVoltage, mediumSweeps[0].sweep.adc0Curr)

#for obj in sensorPackets:
#    print("Acc Z: ", obj.accZ)
#    print('Gyr Z: ', obj.gyroZ)
    
    
    

#for i in range(Burst_Sweep.nSteps):
#    print(round(pckts[1].sweepArray[0].sweepVoltage[i],2), round(pckts[1].sweepArray[0].adc0Curr[i], 2), round(pckts[1].sweepArray[0].adc1Curr[i], 2), round(pckts[1].sweepArray[0].adc2Curr[i], 2))
#for obj in pckts:
#    for i in range(Medium_Sweep.nSteps):
#        print(round(obj.sweep.sweepVoltage[i],2), round(obj.sweep.adc0Curr[i],2), round(obj.sweep.adc1Curr[i],2), round(obj.sweep.adc2Curr[i],2))
#for i in range(pckts[10].sweep.nSteps):
#    print(round(pckts[10].sweep.sweepVoltage[i],2), round(pckts[10].sweep.adcDACVolt[i],2))

#for i in range(pckts[0].nSteps):
#    print(round(pckts[0].sweep.sweepVoltage[i], 2))
     
if __name__ == '__main__':
    #Unit test cases
    packetList = readFile("UDIP0000.DAT")

    #print(packet[0].accX)
    #print(packet[0].accY)
    #print(packet[0].accZ)
    accXList = []
    accYList = []
    accZList = []
    accHList = []
    gyrXList = []
    gyrYList = []
    gyrZList = []
    magXList = []
    magYList = []
    magZList = []
    tmpDList = []
    tmpPList = []
    tmpSList = []
    pd1List = []
    pd2List = []

    sensList = []
    medSwpList = []
    cstSwpList = []
    lrgSwpList = []

    for packet in packetList:
        if(packet.pcktType == typeSens):
            sensList.append(packet)
            accXList.append(packet.accX)
            accYList.append(packet.accY)
            accZList.append(packet.accZ)
            accHList.append(packet.accH)
            gyrXList.append(packet.gyroX)
            gyrYList.append(packet.gyroY)
            gyrZList.append(packet.gyroZ)
            magXList.append(packet.magX)
            magYList.append(packet.magY)
            magZList.append(packet.magZ)
            tmpDList.append(packet.temperature_d)
            tmpPList.append(packet.temperature_p)
            tmpSList.append(packet.temperature_s)
            pd1List.append(packet.pd_1)
            pd2List.append(packet.pd_2)
        elif(packet.pcktType == typeMed):
            medSwpList.append(packet)
        elif(packet.pcktType == typeCst):
            cstSwpList.append(packet)
        elif(packet.pcktType == typeLrg):
            lrgSwpList.append(packet)
        elif(packet.pcktType == typeBrst):
            pass
        else:
            print('Why?')

    #print(medSwpList[0].sweep.vNeg)
    #plt.plot(medSwpList[1].sweep.adcDACVolt,medSwpList[1].sweep.adc0Curr)
    #plt.plot(medSwpList[1].sweep.adcDACVolt,medSwpList[1].sweep.adc1Curr)
    #plt.plot(medSwpList[5].sweep.ProbeVolt,medSwpList[5].sweep.adc2Curr)
    #plt.plot([-6.,6.],[0.,0.], color='k')
    #plt.plot(medSwpList[0].sweep.adcDACVolt)
    #plt.plot(medSwpList[0].sweep.adcDACVal)
    #plt.show()
    #print(medSwpList[0].sweep.adcDACVal)
    #print('\n')
    #print(medSwpList[0].sweep.adc1Curr)

    # plt.plot(lrgSwpList[0].sweep.adcDACVolt,lrgSwpList[0].sweep.adc1Curr)
    # plt.show()
    # print(lrgSwpList[0].sweep.adcDACVal)
    # print('\n')
    # print(lrgSwpList[0].sweep.adc1Curr)

    print(sensList[-1].tFinal - sensList[0].tInitial)

    #print(cstSwpList[6].tFinal - cstSwpList[6].tInitial)
    #print(cstSwpList[0].)
    # print(cstSwpList[0].sweep.ProbeVolt[0])
    # plt.plot(cstSwpList[0].sweep.adc0Curr)
    # plt.plot(cstSwpList[0].sweep.adc1Curr)
    # plt.plot(cstSwpList[0].sweep.adc2Curr)
    # plt.show()

    # print(cstSwpList[1].sweep.ProbeVolt[0])
    # plt.plot(cstSwpList[1].sweep.adc0Curr)
    # plt.plot(cstSwpList[1].sweep.adc1Curr)
    # plt.plot(cstSwpList[1].sweep.adc2Curr)
    # plt.show()

    # print(cstSwpList[2].sweep.ProbeVolt[0])
    # plt.plot(cstSwpList[2].sweep.adc0Curr)
    # plt.plot(cstSwpList[2].sweep.adc1Curr)
    # plt.plot(cstSwpList[2].sweep.adc2Curr)
    # plt.show()

    # print(cstSwpList[3].sweep.ProbeVolt[0])
    # plt.plot(cstSwpList[3].sweep.adc0Curr)
    # plt.plot(cstSwpList[3].sweep.adc1Curr)
    # plt.plot(cstSwpList[3].sweep.adc2Curr)
    # plt.show()

    # print(cstSwpList[4].sweep.ProbeVolt[0])
    # plt.plot(cstSwpList[4].sweep.adc0Curr)
    # plt.plot(cstSwpList[4].sweep.adc1Curr)
    # plt.plot(cstSwpList[4].sweep.adc2Curr)
    # plt.show()

    # print(cstSwpList[5].sweep.ProbeVolt[0])
    # plt.plot(cstSwpList[5].sweep.adc0Curr)
    # plt.plot(cstSwpList[5].sweep.adc1Curr)
    # plt.plot(cstSwpList[5].sweep.adc2Curr)
    # plt.show()

    # print(cstSwpList[6].sweep.ProbeVolt[0])
    # plt.plot(cstSwpList[6].sweep.adc0Curr)
    # plt.plot(cstSwpList[6].sweep.adc1Curr)
    # plt.plot(cstSwpList[6].sweep.adc2Curr)
    # plt.show()

    # print(cstSwpList[7].sweep.ProbeVolt[0])
    # plt.plot(cstSwpList[7].sweep.adc0Curr)
    # plt.plot(cstSwpList[7].sweep.adc1Curr)
    # plt.plot(cstSwpList[7].sweep.adc2Curr)
    # plt.show()

    # print(cstSwpList[8].sweep.ProbeVolt[0])
    # plt.plot(cstSwpList[8].sweep.adc0Curr)
    # plt.plot(cstSwpList[8].sweep.adc1Curr)
    # plt.plot(cstSwpList[8].sweep.adc2Curr)
    # plt.show()

    # print(cstSwpList[9].sweep.ProbeVolt[0])
    # plt.plot(cstSwpList[9].sweep.adc0Curr)
    # plt.plot(cstSwpList[9].sweep.adc1Curr)
    # plt.plot(cstSwpList[9].sweep.adc2Curr)
    # plt.show()
 
    #plt.plot(cstSwpList[3].sweep.adcPD1Volt)
    #plt.plot(cstSwpList[3].sweep.adcDACVolt)
    #print(cstSwpList.sweep.adcDACVolt[1])
    #print([x.sweep.adcDACVal[1] for x in cstSwpList])
    #plt.show()

    

    # Print accelerations
    # print('\n')
    # print('\n')
    # print('Accel X: ',accXList[0:8])
    # print('Accel Y: ',accYList[0:8])
    # print('Accel Z: ',accZList[0:8])
    # print('Accel H: ',accHList[0:8])

    # plt.plot(accXList)
    # plt.plot(accYList)
    # plt.plot(accZList)
    # #plt.plot(accHList)
    # plt.title("Accel")
    # plt.show()

    # # Print gyroscopes
    # print('\n')
    # print('\n')
    # print('Gyro X: ',gyrXList[0:8])
    # print('Gyro Y: ',gyrYList[0:8])
    # print('Gyro Z: ',gyrZList[0:8])

    # plt.plot(gyrXList)
    # plt.plot(gyrYList)
    # plt.plot(gyrZList)
    # plt.title("Gyro")
    # plt.show()

    # # print magnetometer
    # print('\n')
    # print('\n')
    # print('Mag X: ',magXList[0:8])
    # print('Mag Y: ',magYList[0:8])
    # print('Mag Z: ',magZList[0:8])

    # plt.plot(magXList)
    # plt.plot(magYList)
    # plt.plot(magZList)
    # plt.title("Mag")
    # plt.show()
    #print(magYList[20:30])
    #print([unpack('<h', x.pyld[16:18])[0] for x in sensList[20:30]])
    #unpack('<h', sensList.pyld[16:18])
    # print('Mag X: ',magXList)
    # print('Mag Y: ',magYList)
    # print('Mag Z: ',magZList)

    # # print temps
    # print('\n')
    # print('\n')
    # print('Temp digital: ',tmpDList[0:8])
    # print('Temp port: ',tmpPList[0:8])
    # print('Temp skin: ',tmpSList[0:8])

    # plt.plot(tmpDList)
    # plt.plot(tmpPList)
    # plt.plot(tmpSList)
    # plt.title("Temp")
    # plt.show()

    # # print photodiodes
    # print('\n')
    # print('\n')
    # print('PD1: ',pd1List[0:8])
    # print('PD2: ',pd2List[0:8])

    # plt.plot(pd1List)
    # plt.plot(pd2List)
    # plt.title("PD")
    # plt.show()
    #print(pd1List[0:100])