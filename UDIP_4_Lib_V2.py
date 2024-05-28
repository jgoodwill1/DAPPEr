# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 20:35:53 2020
Class definitions for Reading Packets from University of Delaware Ionosphere
Probe RockSat-C 2020 Mission
Last updated 2024/05/20
@author: Jarod Dagney
        Jarrod Bieber
"""
from numpy import tile

from struct import unpack

from scipy import interpolate

#import matplotlib.pyplot as plt

#import matplotlib.pyplot as plt
#from scipy.io import savemat
#import numpy as np

VERBOSE = False
lenHedr = 15

lenSens = 24

lenBrst = 133 * 10 * 2 * 4 + 4
lenMed = 265 * 2 * 4 + 4
lenLrg = 253 * 2 * 4 + 4
lenFull = 256 * 2 * 4 + 8
lenDense = 256 * 2 * 4 + 8

typeSens = 0x01
typeMed  = 0x09 # 0x10
typeLrg  = 0x11
typeBrst = 0x20
typeFull_Probe = 0x30  # Full/Linear Sweep - 0x40
typeFull_Rckt = 0x50  # Full/Linear Sweep - Rocket Ground
typeDens_Probe = 0x10  # High Density/Special Sweep - 0x60
typeDens_Rckt = 0x70  # High Density/Special Sweep - Rocket Ground

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
    Large_Sweep, Burst_Sweep, Full_Sweep, Dense_Sweep.
    
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
        accScale: Maximum range of accelerometer 2, 4, 8, or 16 g
        accX, accY, accZ: Acceleration in g (9.81 m/s^2) in three axes, from 9 Deg. of Freedom Sensor
        gyrX, gyrY, gyrZ: Spin rate in degrees per second in three axes
        magX, magY, magZ: Magnetic field in Gauss in three axes
        accHiRng: High range accelerometer, Z-direction (from accel)
        temperature: Internal canister temperature
        phDiode: Photo Diode Voltage (?)
    """
    SENSITIVITY_ACCELEROMETER_2 = 0.000061
    SENSITIVITY_ACCELEROMETER_4 = 0.000122
    SENSITIVITY_ACCELEROMETER_8 = 0.000244
    SENSITIVITY_ACCELEROMETER_16= 0.000732
    
    SENSITIVITY_GYROSCOPE_2000  = 0.07
    
    SENSITIVITY_MAGNETOMETER_4  = 0.00014
    
    temperatureConversion = 1/8.#0.0625

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
        self.accScale = unpack('<B', self.pyld[lenSens - 1:lenSens])[0]
        #self.acc = unpack('<hhhH',self.pyld[0:8])#three signed two byte integers and 1 unsigned two byte integer
        self.accMed = unpack('<hhh',self.pyld[0:6])
        self.accHiRng = unpack('<H',self.pyld[6:8])  # Do something with this?
        self.calibrateDigAcc()
        #Acceleration data from analog sensor
        #self.accAna =(self.acc[3] * (3.3/4095.) * 2 - 1.30) / (0.02)
        
        #Unpack Gyroscope data from payload
        self.gyr = unpack('<hhh', self.pyld[8:14])
        self.calibrateDigGyr()
        
        #Unpack Magnetic Field data from payload
        self.mag = unpack('<hhh', self.pyld[14:20])
        self.calibrateDigMag()
        
        #Unpack Temperature data from payload
        tmpD = unpack('<H', self.pyld[20:22])[0]
        self.phDiode = unpack('<H', self.pyld[22:24])[0]  # Check if this is okay
        self.phDiode = self.phDiode * 3.3 / 4095 / 16
        #self.temperature = tmpD* self.temperatureConversion
        self.temperature = ((tmpD * 3300 / 4095) - 500) / 10
        #self.temperatureAna = ((tmpA * (3.3/4095.))-0.5)*100
        
            
        
    def calibrateDigAcc(self):
        """Calibrates the acceleration and stores it in respective class attributes in g"""
        # if(self.accScale == 2):
        #     self.accX = self.accMed[0] * self.SENSITIVITY_ACCELEROMETER_2
        #     self.accY = self.accMed[1] * self.SENSITIVITY_ACCELEROMETER_2
        #     self.accZ = self.accMed[2] * self.SENSITIVITY_ACCELEROMETER_2
        # elif(self.accScale == 4):
        #     self.accX = self.accMed[0] * self.SENSITIVITY_ACCELEROMETER_4
        #     self.accY = self.accMed[1] * self.SENSITIVITY_ACCELEROMETER_4
        #     self.accZ = self.accMed[2] * self.SENSITIVITY_ACCELEROMETER_4
        # elif(self.accScale == 8):
        #     self.accX = self.accMed[0] * self.SENSITIVITY_ACCELEROMETER_8
        #     self.accY = self.accMed[1] * self.SENSITIVITY_ACCELEROMETER_8
        #     self.accZ = self.accMed[2] * self.SENSITIVITY_ACCELEROMETER_8
        # elif(self.accScale == 16):
        #     self.accX = self.accMed[0] * self.SENSITIVITY_ACCELEROMETER_16
        #     self.accY = self.accMed[1] * self.SENSITIVITY_ACCELEROMETER_16
        #     self.accZ = self.accMed[2] * self.SENSITIVITY_ACCELEROMETER_16
        # else:
        #     #Not valid accel scale
        #     self.accX = 0;
        #     self.accY = 0;
        #     self.accZ = 0;

        # self.accX = self.accMed[0] * self.SENSITIVITY_ACCELEROMETER_16
        # self.accY = self.accMed[1] * self.SENSITIVITY_ACCELEROMETER_16
        # self.accZ = self.accMed[2] * self.SENSITIVITY_ACCELEROMETER_16
        self.accX = self.accMed[0] / 95.43
        self.accY = self.accMed[1] / 95.43
        self.accZ = self.accMed[2] / 95.43
        
    
    def calibrateDigGyr(self):
        """Calibrates the spin rate and stores it in respective class attributes in degrees per second"""
        # self.gyroX = self.gyr[0] * self.SENSITIVITY_GYROSCOPE_2000
        # self.gyroY = self.gyr[1] * self.SENSITIVITY_GYROSCOPE_2000
        # self.gyroZ = self.gyr[2] * self.SENSITIVITY_GYROSCOPE_2000
        self.gyroX = self.gyr[0] / 936.25
        self.gyroY = self.gyr[1] / 936.25
        self.gyroZ = self.gyr[2] / 936.25
        
    
    def calibrateDigMag(self):
        """Calibrates the magnetic field data and stores it in respective class attributes in Gauss"""
        # self.magX = self.mag[0] * self.SENSITIVITY_MAGNETOMETER_4
        # self.magY = self.mag[1] * self.SENSITIVITY_MAGNETOMETER_4
        # self.magZ = self.mag[2] * self.SENSITIVITY_MAGNETOMETER_4
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


class Full_Sweep(Packet):
    """
        Subclass of Packet, creates an object of type Sweep.

        Attributes:
            All Packet Attributes
            sweep: Object of type Sweep that holds Voltage and Current values
            grndType: '0' for grounded to rocket skin, '1' for disconnected from rocket skin
     """
    nSteps = 256

    def __init__(self, count, tInitial, tFinal, pcktType, payloadLen, grndType, pyld):
        super().__init__(count, tInitial, tFinal, pcktType, payloadLen)
        self.sweep = Sweep(pyld, self.nSteps)
        self.grndType = grndType
        return


class Dense_Sweep(Packet):
    """
        Subclass of Packet, creates an object of type Sweep.

        Attributes:
            All Packet Attributes
            sweep: Object of type Sweep that holds Voltage and Current values
            grndType: '0' for grounded to rocket skin, '1' for disconnected from rocket skin
     """
    nSteps = 256

    def __init__(self, count, tInitial, tFinal, pcktType, payloadLen, grndType, pyld):
        super().__init__(count, tInitial, tFinal, pcktType, payloadLen)
        self.sweep = Sweep(pyld, self.nSteps)
        self.grndType = grndType
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
    nStepsFull = 256  # New for UDIP 4
    nStepsHiDens = 256  # New for UDIP 4
    
    nAvg = 4.0  # 2.0? 10.0? 4.0? 1.0?

    # No Longer Used to get Current in Transfer Functions
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
    #RP1 = 49900.0
    RP1 = 20000.0
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

        # Get new equation from Ranyah
        # Sample 16 times(?), then formula
        #self.vPos = (self.vPosVal / 16.0) * (3.3 / 4096.0) * (self.RP1 + self.RP2) / self.RP2
        #self.vPos = (self.vPosVal / 3.0)
        #self.vPos = (self.vPosVal / self.nAvg) * (4.75 / 1023) * (self.RP1 + self.RP2) / self.RP2
        self.vPos = (self.vPosVal / self.nAvg) * (4.75 / 1023) * (self.RP1 + self.RP2) / self.RP2
        #self.vPos = (self.vPosVal / self.nAvg) * (5.0 / 1023) * (self.RP1 + self.RP2) / self.RP2
        # Get new equation from Ranyah
        #self.vNeg = (self.RN1 + self.RN2) / self.RN1 * ( ( (self.vNegVal / 16.0) * (3.3 / 4096.0) ) - (self.vPos ) * (self.RN2 / (self.RN1 + self.RN2)))
        self.vNeg = - (self.vNegVal / 3.0)

        self.phDiodeInitial = unpack('<H', pyld[4:6])[0]
        self.phDiodeInitial = self.phDiodeInitial * 3.3 / 4095 / 16

        self.phDiodeFinal = unpack('<H', pyld[6:8])[0]
        self.phDiodeFinal = self.phDiodeFinal * 3.3 / 4095 / 16

        
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
            arr = unpack('<HHHH', self.pyld[(4*2*i) + 4 + 4: (4*2*i + 8) + 4 + 4])
            self.adcDACVal[i] = arr[0] / self.nAvg
            #print(arr[0])
            self.adc0Val[i] = arr[1] / self.nAvg
            #print(arr[1])
            self.adc1Val[i] = arr[2] / self.nAvg
            self.adc2Val[i] = arr[3] / self.nAvg
            
        return

    
    def fillVolt(self):
        """Converts a value from ADC to the voltage"""
        conv = 4.75 / 1023  # 5.0 / 4095.  # need new conversion factor - 5.0 / 1023, 4.76 / 1023, 4.0 / 1023
        for i in range(self.nSteps):
            self.adcDACVolt[i] = self.adcDACVal[i] * conv
            self.adc0Volt[i] = self.adc0Val[i] * conv
            self.adc1Volt[i] = self.adc1Val[i] * conv
            self.adc2Volt[i] = self.adc2Val[i] * conv
        return
         
    def fillCurr(self):
        """"Converts the voltage to ADC to a current from the probe"""
        for i in range(self.nSteps):
           #self.adc0Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc0Volt[i] - (self.R1/(self.R1+self.R2)) * 11.97 ) * ((self.RG1/(self.RG1+self.RG2))**0) * 1E9
           #self.adc1Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc1Volt[i] - (self.R1/(self.R1+self.R2)) * 11.97 ) * ((self.RG1/(self.RG1+self.RG2))**1) * 1E9
           #self.adc2Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc2Volt[i] - (self.R1/(self.R1+self.R2)) * 11.97 ) * ((self.RG1/(self.RG1+self.RG2))**2) * 1E9

            #self.adc0Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc0Volt[i] - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**0) * 1E9
            #self.adc1Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc1Volt[i] - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**1) * 1E9
            #self.adc2Curr[i] = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * self.adc2Volt[i] - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**2) * 1E9

           #self.adc0Curr[i] = (-6 * (self.adc0Volt[i] - 2) / 5e5
           #self.adc0Curr[i] = (-6 * self.adc0Volt[i] + 12) / 500000
           self.adc0Curr[i] = (-6 * self.adc0Volt[i] + self.vPos) / 500000

           #self.adc1Curr[i] = (-6 * (self.adc1Volt[i] - 2) / 1e7)
           #self.adc1Curr[i] = (-6 * self.adc1Volt[i] + 12) / 1e7
           self.adc1Curr[i] = (-6 * self.adc1Volt[i] + self.vPos) / 1e7

           #self.adc2Curr[i] = (-6 * (self.adc2Volt[i] - 2) / 2e8)
           #self.adc2Curr[i] = (-6 * self.adc2Volt[i] + 12) / 2e8
           self.adc2Curr[i] = (-6 * self.adc2Volt[i] + self.vPos) / 2e8

        #self.adc0CurrMax = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 3.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**0) * 1E9
        #self.adc1CurrMax = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 3.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**1) * 1E9
        #self.adc2CurrMax = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 3.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**2) * 1E9
        #self.adc0CurrMin = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 0.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**0) * 1E9
        #self.adc1CurrMin = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 0.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**1) * 1E9
        #self.adc2CurrMin = (1/self.RF) * ((self.R1+self.R2)/self.R2) * ( ((self.R4+self.R3)/self.R4) * 0.0 - (self.R1/(self.R1+self.R2)) * self.vPos ) * ((self.RG1/(self.RG1+self.RG2))**2) * 1E9
        self.adc0CurrMax = (-6 * (4 - 2) / 5e5)
        self.adc1CurrMax = (-6 * (4 - 2) / 1e7)
        self.adc2CurrMax = (-6 * (4 - 2) / 2e8)
        self.adc0CurrMin = (-6 * (0 - 2) / 5e5)
        self.adc1CurrMin = (-6 * (0 - 2) / 1e7)
        self.adc2CurrMin = (-6 * (0 - 2) / 2e8)

        '''
        ADC 0 Min = -24937.5 nA, Max = 22062.5 nA
        ADC 1 Min = -1234.1 nA, Max = 1050.6 nA
        ADC 2 Min = -58.8 nA, Max = 50.03 nA
        '''

        return
           
        
    def sweepTransfer(self):
        """Converts voltage to ADC from DAC pin to voltage applied to probe"""
        #self.sweepVoltage = (self.RS3/self.RS1 + self.RS3/self.RS2 + 1) * self.adcSweepVolt - (self.RS3/self.RS1) * 11.97
        #self.sweepVoltage = ( ( self.RS1 + self.RS2 ) / self.RS2 * self.adcSweepVolt) - ( self.RS1 / self.RS2 ) * self.vNeg
        #self.sweepVoltage = (self.RS3/self.RS1 + self.RS3/self.RS2 + 1) * self.adcDACVolt - (self.RS3/self.RS1) * self.vPos
        self.sweepVoltage = (self.adcDACVolt * 6) - 12  # Check if this is correct
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
    #raw = [:
    loc = 0
    i = 0
    while(loc < len(raw)):
        #print()
        #print(raw[loc:loc+15])
        #print()
        sync     = unpack('<BB',raw[loc + 0:loc + 2])
        count    = unpack('<H', raw[loc + 2:loc + 4])[0]
        tInitial = unpack('<I', raw[loc + 4:loc + 8])[0]
        tFinal   = unpack('<I', raw[loc + 8:loc + 12])[0]
        pcktType = unpack('<B', raw[loc + 12:loc + 13])[0]
        pyldLen  = unpack('<H', raw[loc + 13:loc + 15])[0]

        #print(f'sync = {sync}')
        #print(f'count = {count}')
        #print(f'tInitial = {tInitial}')
        #print(f'tFinal = {tFinal}')
        #print(f'pcktType = {pcktType}')
        #print(f'pyldLen = {pyldLen}')
        #print(raw[loc+pyldLen-20:loc+pyldLen+20])
        #print(raw[loc+pyldLen:loc+pyldLen+1])

        if(verifyHeader(sync, pcktType, pyldLen)):
            print("valid header YAY")
            if(pcktType == typeSens):
                packet = Sensor_Packet(count, tInitial, tFinal, pcktType, pyldLen, raw[(loc+lenHedr) : (loc+lenHedr+pyldLen)] )
                myPackets.append(packet)
            elif(pcktType == typeFull_Rckt):
                print(f'Packet type Full - Rocket')
                print(f'Packet {i}, loc = {loc}')
                packet = Full_Sweep(count, tInitial, tFinal, pcktType, pyldLen, 0, raw[(loc+lenHedr) : (loc+lenHedr+pyldLen)])
                myPackets.append(packet)
            elif(pcktType == typeFull_Probe):
                print(f'Packet type Full')
                print(f'Packet {i}, loc = {loc}')
                packet = Full_Sweep(count, tInitial, tFinal, pcktType, pyldLen, 1, raw[(loc + lenHedr): (loc + lenHedr + pyldLen)])
                myPackets.append(packet)
            elif (pcktType == typeDens_Rckt):
                print(f'Packet type Dense - Rocket')
                print(f'Packet {i}, loc = {loc}')
                packet = Dense_Sweep(count, tInitial, tFinal, pcktType, pyldLen, 0, raw[(loc + lenHedr): (loc + lenHedr + pyldLen)])
                myPackets.append(packet)
            elif (pcktType == typeDens_Probe):
                print(f'Packet type Dense')
                print(f'Packet {i}, loc = {loc}')
                packet = Dense_Sweep(count, tInitial, tFinal, pcktType, pyldLen, 1, raw[(loc + lenHedr): (loc + lenHedr + pyldLen)])
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
            else:
                loc+=1
            loc += lenHedr+pyldLen
            if(VERBOSE):
                print("Packet #: ", i)
            i+=1
        else:
            loc+=1
    print()
    print(f'Generated {i} packets')
    print()
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
        print(f'sync[0] = {sync[0]} (should be 0x55), sync[1] = {sync[1]} (should be 0x44)')
        return False
    if(pcktType == typeSens):
        print("packet type Sensor")
        if(pyldLen != lenSens):
            print(f"bad length: {pyldLen} instead of {lenSens}")
            return False
        else:
            print("good length")
            return True
    if (pcktType == typeFull_Probe):
        print("packet type Full")
        if (pyldLen != lenFull):
            print(f"bad length: {pyldLen} instead of {lenFull}")
            return False
        else:
            print("good length")
            return True
    if (pcktType == typeFull_Rckt):
        print("packet type Full - Rocket")
        if (pyldLen != lenFull):
            print(f"bad length: {pyldLen} instead of {lenFull}")
            return False
        else:
            print("good length")
            return True
    if (pcktType == typeDens_Probe):
        print("packet type Dense")
        if (pyldLen != lenDense):
            print(f"bad length: {pyldLen} instead of {lenDense}")
            return False
        else:
            print("good length")
            return True
    if (pcktType == typeDens_Rckt):
        print("packet type Dense - Rocket")
        if (pyldLen != lenDense):
            print(f"bad length: {pyldLen} instead of {lenDense}")
            return False
        else:
            print("good length")
            return True
    if(pcktType == typeMed):
        print("packet type Med???")
        if(pyldLen != lenMed):
            return False
        else:
            return True
    if(pcktType == typeLrg):
        print("packet type Large???")
        if(pyldLen != lenLrg):
            return False
        else:
            return True
    if(pcktType == typeBrst):
        print("packet type Burst???")
        if(pyldLen != lenBrst):
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


def findIndexs(mypackets): #find indexs of various packet types
    sensor = []
    full = []
    dense = []
    #medium = []
    #large = []
    #burst = []
    i = 0
    #pcktType
    typeSens = 0x01
    #typeMed = 0x10
    #typeLrg = 0x11
    #typeBrst = 0x20
    typeFull_Probe = 0x30  # Full/Linear Sweep - 0x40
    typeFull_Rckt = 0x50  # Full/Linear Sweep - Rocket Ground
    typeDens_Probe = 0x10  # High Density/Special Sweep - 0x60
    typeDens_Rckt = 0x70  # High Density/Special Sweep - Rocket Ground

    for obj in mypackets:
        if(obj.pcktType == typeSens):
            sensor.append(i)
        elif (obj.pcktType == typeFull_Probe):
            full.append(i)
        elif (obj.pcktType == typeFull_Rckt):
            full.append(i)
        elif (obj.pcktType == typeDens_Probe):
            dense.append(i)
        elif (obj.pcktType == typeDens_Rckt):
            dense.append(i)
        #elif(obj.pcktType == typeMed):
        #    medium.append(i)
        #elif(obj.pcktType == typeLrg):
        #    large.append(i)
        #elif(obj.pcktType == typeBrst):
        #    burst.append(i)
        i=i+1
    #return sensor, medium, large, burst
    return sensor, full, dense
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

# pckts = readFile('UDIP0013.DAT')

# temps = []

# for p in pckts:
#     if p.pcktType == typeSens:
#         temps.append(p.temperature)
        
# plt.plot(temps)


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


packets = readFile('UDIP100.dat')
sensorIndex, fullIndex, denseIndex = findIndexs(packets)

sweep1 = packets[fullIndex[3]].sweep
dacVal = sweep1.adcDACVal
swpVolt = sweep1.sweepVoltage

sweep2 = packets[denseIndex[0]].sweep
dacVal_2 = sweep2.adcDACVal
swpVolt2 = sweep2.sweepVoltage

print(f'dacVal (full) = {dacVal}')
print()
print(f'dacVal (dense) = {dacVal_2}')

print(f'sweepVoltage (full) = {swpVolt}')
print()
print(f'sweepVoltage (dense) = {swpVolt2}')

# for i in range(len(dacVal)):
#     if(dacVal[i] == dacVal_2[i]):
#         print('God has forsaken us!!!')
