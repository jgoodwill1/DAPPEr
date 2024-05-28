"""
Created on November 29, 2023

Last Updated on May 13th, 2024

@author: Jarod Dagney
        Jarrod Bieber
"""

#import numpy as np
#import os
#from struct import unpack

#from scipy import signal

#import pickle

#import matplotlib.pyplot as plt

#import UDIP_Lib_V19 as UDIP_Lib
import UDIP_4_Lib_V2 as UDIP_Lib  #making updating UDIP_Lib easier

#import RockSat_1_14_fitting_functions as R_fitting

import json

#from datetime import datetime


#global variables

def export_data():
    '''
    create a json file with keys  packet : [ Sensor, Full, Dense ]
    Each Full, Dense keys have values of their index ranging from 0 to the number of that type of packet
    For each index, there are keys:
        Start - start time of this packet after launch in ms
        Stop - stop time of this packet after launch in ms
        Count - the index of this packet among all packets
        GroundType - '0' for grounded to rocket skin, '1' for disconnected from rocket skin
        Voltage - a list of sweep voltages
        Current0 - a list of measured currents from the 0 gain stage
        Current1 - a list of measured currents from the 1 gain stage
        Current2 - a list of measured currents from the 2 gain stage
    '''
    #mypackets = UDIP_Lib.readFile("UDIP0016.DAT")
    mypackets = UDIP_Lib.readFile("UDIP115.DAT")  # Update this file name
    #ind_sensor, ind_med, ind_large, ind_burst = findIndexs(mypackets)
    ind_sensor, ind_full, ind_dense = findIndexs(mypackets)

    #json_obj = { "packets" : {"Sensor" : {}, "Medium" : {}, "Large" : {}, "Burst" : {}} }
    json_obj = {"packets": {"Sensor": {}, "Full": {}, "Dense": {}}}
    print("Adding sensor packets")
    print(f'Index Sensor Length = {len(ind_sensor)}')
    for i in range(len(ind_sensor)):
        pckt = mypackets[ind_sensor[i]]
        json_obj["packets"]["Sensor"][i] = {'Start': pckt.tInitial, 'Stop': pckt.tFinal, 'Count': pckt.count,
                                            #"Acceleration": [pckt.accX, pckt.accY, pckt.accZ, pckt.accAna],
                                            "Acceleration": [pckt.accX, pckt.accY, pckt.accZ, pckt.accHiRng],
                                            "Spin Rate": [pckt.gyroX, pckt.gyroY, pckt.gyroZ],
                                            "Magnetic Field": [pckt.magX, pckt.magY, pckt.magZ],
                                            #"Temperature": [pckt.temperature, pckt.temperatureAna]}
                                            "Temperature": [pckt.temperature],
                                            "Photo Diode": [pckt.phDiode]}
    print("Adding full sweep packets")
    print(f'Index Full Length = {len(ind_full)}')
    for i in range(len(ind_full)):
        pckt = mypackets[ind_full[i]]
        #print()
        #print(f'vPosVal = {pckt.sweep.vPosVal}')
        #print(f'vNegVal = {pckt.sweep.vNegVal}')
        # print(f'vPos = {pckt.sweep.vPos}')
        # print(f'vNeg = {pckt.sweep.vNeg}')
        # print(f'adcDACVal = {pckt.sweep.adcDACVal}')
        # print(f'adcDACVolt = {pckt.sweep.adcDACVolt}')
        # print(f'sweepVoltage = {pckt.sweep.sweepVoltage}')
        # print(f'adc0Val = {pckt.sweep.adc0Val}')
        # print(f'adc0Volt = {pckt.sweep.adc0Volt}')
        # print(f'adc0Curr = {pckt.sweep.adc0Curr}')
        # print(f'adc1Val = {pckt.sweep.adc1Val}')
        # print(f'adc1Volt = {pckt.sweep.adc1Volt}')
        # print(f'adc1Curr = {pckt.sweep.adc1Curr}')
        # print(f'adc2Val = {pckt.sweep.adc2Val}')
        # print(f'adc2Volt = {pckt.sweep.adc2Volt}')
        # print(f'adc2Curr = {pckt.sweep.adc2Curr}')
        # print(f'adc0CurrMax = {pckt.sweep.adc0CurrMax}')
        # print(f'adc1CurrMax = {pckt.sweep.adc1CurrMax}')
        # print(f'adc2CurrMax = {pckt.sweep.adc2CurrMax}')
        # print(f'adc0CurrMin = {pckt.sweep.adc0CurrMin}')
        # print(f'adc1CurrMin = {pckt.sweep.adc1CurrMin}')
        # print(f'adc2CurrMin = {pckt.sweep.adc2CurrMin}')
        # print()
        json_obj["packets"]["Full"][i] = {'Start': pckt.tInitial, 'Stop': pckt.tFinal, 'Count': pckt.count,
                                          'GroundType': pckt.grndType, 'Voltage': pckt.sweep.sweepVoltage.tolist(),
                                          'Current0': pckt.sweep.adc0Curr.tolist(),
                                          'Current1': pckt.sweep.adc1Curr.tolist(),
                                          'Current2': pckt.sweep.adc2Curr.tolist(),
                                          'PhotoDiodeInitial': pckt.sweep.phDiodeInitial,
                                          'PhotoDiodeFinal': pckt.sweep.phDiodeFinal}
    print("Adding dense sweep packets")
    print(f'Index Dense Length = {len(ind_dense)}')
    for i in range(len(ind_dense)):
        pckt = mypackets[ind_dense[i]]
        json_obj["packets"]["Dense"][i] = {'Start': pckt.tInitial, 'Stop': pckt.tFinal, 'Count': pckt.count,
                                          'GroundType': pckt.grndType, 'Voltage': pckt.sweep.sweepVoltage.tolist(),
                                          'Current0': pckt.sweep.adc0Curr.tolist(),
                                          'Current1': pckt.sweep.adc1Curr.tolist(),
                                          'Current2': pckt.sweep.adc2Curr.tolist(),
                                          'PhotoDiodeInitial': pckt.sweep.phDiodeInitial,
                                          'PhotoDiodeFinal': pckt.sweep.phDiodeFinal}

#    print("Adding medium sweep packets")
#    for i in range(len(ind_med)):
#        pckt = mypackets[ind_med[i]]
#        json_obj["packets"]["Medium"][i] = {'Start' : pckt.tInitial, 'Stop' : pckt.tFinal, 'Count' : pckt.count, 'Voltage' : pckt.sweep.sweepVoltage.tolist(),
#                                            'Current0' : pckt.sweep.adc0Curr.tolist(), 'Current1' : pckt.sweep.adc1Curr.tolist(), 'Current2' : pckt.sweep.adc2Curr.tolist() }
#    print("Adding large sweep packets")
#    for i in range(len(ind_large)):
#        pckt = mypackets[ind_large[i]]
#        json_obj["packets"]["Large"][i] = {'Start' : pckt.tInitial, 'Stop' : pckt.tFinal, 'Count' : pckt.count, 'Voltage' : pckt.sweep.sweepVoltage.tolist(),
#                                            'Current0' : pckt.sweep.adc0Curr.tolist(), 'Current1' : pckt.sweep.adc1Curr.tolist(), 'Current2' : pckt.sweep.adc2Curr.tolist() }
#    print("Adding burst sweep packets")
#    for i in range(len(ind_burst)):
#        pckt = mypackets[ind_burst[i]]
#        json_obj['packets']['Burst'][i] = {}
#        for j in range(10):
#            json_obj["packets"]["Burst"][i][j] = {'Start' : pckt.tInitial, 'Stop' : pckt.tFinal, 'Count' : pckt.count, 'Voltage' : pckt.sweepArray[j].sweepVoltage.tolist(),
#                                                'Current0' : pckt.sweepArray[j].adc0Curr.tolist(), 'Current1' : pckt.sweepArray[j].adc1Curr.tolist(), 'Current2' : pckt.sweepArray[j].adc2Curr.tolist() }
    json_fileName = "2024_packet_FMS_1.json"

    json_file = open(json_fileName, 'w')
    json.dump(json_obj, json_file, indent=1)
    json_file.close()
    print(f"Wrote to {json_fileName}")
    return


def import_data(fileName):
    '''
    Function to read in json formatted date created from export_data()
    '''
    json_file = open(fileName, 'r')
    json_obj = json.load(json_file)
    json_file.close()

    return json_obj


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
        if(obj.pcktType != typeSens):
            print(f'Packet Type = {obj.pcktType}')
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


#How to parse data
#   Assuming you called packets = import_data(2021_packet.json)
#   1. to get the first medium sweep 
#       first_sweep = packets['packets']['Medium']['0']
#
#   first_sweep containts "Start", "Stop", "Count", "Voltage", "Current0", "Current1", "Current2"
#   2. to then plot the current versus voltage, access the lists in the "Voltage" and "Current#"
#   sweep_voltage = first_sweep["Voltage"]
#   current0 = first_sweep["Current0"]
#   current1 = first_sweep["Current1"]
#   current2 = first_sweep["Current2"]
#       
#

# packets = import_data('../2022/2022_packet.json')
# first_sweep = packets['packets']['Medium']['0']
# print(first_sweep)

export_data()
