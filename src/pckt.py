import numpy as np
import os
import json
import pandas as pd
import matplotlib.pyplot as plt

class pck:
    def __init__(self, path, num, sw_type):
        self.path = path
        self.num = num
        self.sw_type = sw_type
        data = pd.read_json(path)
        data = data['packets'][sw_type][f'{self.num}']
        self.data = data
        self.start = data['Start']/1000
        self.stop  = data['Stop']/1000
        self.grnd = data['GroundType']
        self.iv = pd.DataFrame({
                    'V': np.array(self.data['Voltage'][6:]), 
                    'I0': -np.array(self.data['Current0'][6:]), 
                    'I1': -np.array(self.data['Current1'][6:]), 
                    'I2': -np.array(self.data['Current2'][6:])
                    })
    def VP(self):
        self.iv['dI_dV'] = np.abs(self.iv['I0'].diff())
        ind = self.iv['dI_dV'].argmax()
        return self.iv['V'][ind]
    def Vf(self):
        sign = np.where((np.diff(np.sign(self.iv['I0'])) !=0))[0]
        low = np.argmin(np.abs(self.iv['V'][sign]))
        Vf = self.iv['V'][sign[low]]
        return Vf