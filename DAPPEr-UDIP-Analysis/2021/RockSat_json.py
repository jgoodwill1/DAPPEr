"""
Created on November 29, 2023

@author: Jarod Dagney
"""

import numpy as np
import os
#from struct import unpack

from scipy import signal

#import pickle

#import matplotlib.pyplot as plt

import UDIP_Lib_V19 as UDIP_Lib  #making updating UDIP_Lib easier

import RockSat_1_14_fitting_functions as R_fitting

import json

#from datetime import datetime


#global variables

def export_data():
    '''
    create a json file with keys  packet : [ Sensor, Medium, Large, Burst ]
    Each Medium, Large, Burst keys have values of their index ranging from 0 to the number of that type of packet
    Burst sweep has an addition key for the index of the sweep ranging from 0-10
    For each index, there are keys:
        Start - start time of this packet after launch in ms
        Stop - stop time of this packet after launch in ms
        Count - the index of this packet among all packets
        Voltage - a list of sweep voltages
        Current0 - a list of measured currents from the 0 gain stage
        Current1 - a list of measured currents from the 1 gain stage
        Current2 - a list of measured currents from the 2 gain stage
    '''
    mypackets = UDIP_Lib.readFile("UDIP0016.DAT")
    ind_sensor, ind_med, ind_large, ind_burst = R_fitting.findIndexs(mypackets)

    json_obj = { "packets" : {"Sensor" : {}, "Medium" : {}, "Large" : {}, "Burst" : {}} }
    print("Adding sensor packets")
    for i in range(len(ind_sensor)):
        pckt = mypackets[ind_sensor[i]]
        json_obj["packets"]["Sensor"][i] = {'Start' : pckt.tInitial, 'Stop' : pckt.tFinal, 'Count' : pckt.count, "Acceleration" : [pckt.accX, pckt.accY, pckt.accZ, pckt.accAna],
                                            "Spin Rate" : [pckt.gyroX, pckt.gyroY, pckt.gyroZ], "Magnetic Field" : [pckt.magX, pckt.magY, pckt.magZ], "Temperature" : [pckt.temperature, pckt.temperatureAna]}
    print("Adding medium sweep packets")
    for i in range(len(ind_med)):
        pckt = mypackets[ind_med[i]]
        json_obj["packets"]["Medium"][i] = {'Start' : pckt.tInitial, 'Stop' : pckt.tFinal, 'Count' : pckt.count, 'Voltage' : pckt.sweep.sweepVoltage.tolist(),
                                            'Current0' : pckt.sweep.adc0Curr.tolist(), 'Current1' : pckt.sweep.adc1Curr.tolist(), 'Current2' : pckt.sweep.adc2Curr.tolist() }
    print("Adding large sweep packets")
    for i in range(len(ind_large)):
        pckt = mypackets[ind_large[i]]
        json_obj["packets"]["Large"][i] = {'Start' : pckt.tInitial, 'Stop' : pckt.tFinal, 'Count' : pckt.count, 'Voltage' : pckt.sweep.sweepVoltage.tolist(),
                                            'Current0' : pckt.sweep.adc0Curr.tolist(), 'Current1' : pckt.sweep.adc1Curr.tolist(), 'Current2' : pckt.sweep.adc2Curr.tolist() }
    print("Adding burst sweep packets")
    for i in range(len(ind_burst)):
        pckt = mypackets[ind_burst[i]]
        json_obj['packets']['Burst'][i] = {}
        for j in range(10):
            json_obj["packets"]["Burst"][i][j] = {'Start' : pckt.tInitial, 'Stop' : pckt.tFinal, 'Count' : pckt.count, 'Voltage' : pckt.sweepArray[j].sweepVoltage.tolist(),
                                                'Current0' : pckt.sweepArray[j].adc0Curr.tolist(), 'Current1' : pckt.sweepArray[j].adc1Curr.tolist(), 'Current2' : pckt.sweepArray[j].adc2Curr.tolist() }
    json_fileName = "2021_packet.json"
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
