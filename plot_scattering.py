# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:02:46 2019

Comparison between the Jason-2 sigma0 and the NRCS from the OceanSAR simulator

@author: lyh
"""

# plot scattering cahnges as a function of wind speed

# Setting up the packages
import pickle
import numpy as np
from matplotlib import pyplot as plt
import os


wind_speed=[1,2,3,5,7,9,11,13,15,20,25]
x = np.linspace(1,25,25)
data = np.zeros((np.size(wind_speed), 64))
index = np.zeros((np.size(wind_speed), 64))
file_address = r"D:\research\TU Delft\programs\Alticube_tool\Alticubes'tool\scattering"



for i in range(np.size(wind_speed)):
    name = 'NRCS' + '_' + str(wind_speed[i]) + '.txt'
    path = file_address + os.sep + name
    f=open(path, 'rb')
    d=pickle.load(f)
    data[i,:] = 10 * np.log10(d)
    index[i,:] = np.ones_like(d) * wind_speed[i]
data = data.reshape((data.size,))
index = index.reshape((index.size,))
plt.figure()
f = np.polyfit(index, data, 9)
fit = np.poly1d(f)
ll = plt.plot(wind_speed, fit(wind_speed), label="Fitting curve")
lo = plt.scatter(index, data,marker='x', color=(0.8,0.,0.), label="Samples") 
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('Wind speed (m/s)', fontsize=16)
plt.ylabel('NRCS (dB)', fontsize=16)
plt.legend(loc="best", fontsize=16)
plt.grid()