# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 12:11:32 2019

It is the program to calculate the height ambiguity in the concept 3

@author: lyh
"""
# Setting up the packages
import numpy as np
from matplotlib import pyplot as plt
from constants import constants as const
from orbit.orbit import orbit_com
import ocs_io as tpio

# Setting up the paths
cfg_file = r"D:\research\TU Delft\programs\ALTCUBTOOL\cfg\parameters.cfg"
cfg = tpio.ConfigFile(cfg_file)
baseline = np.linspace(0.5,1000, 1000)

# lambda [m]
l0 = const.c / cfg.radar.f0
# along tr. bw [rad]
a_bw = l0 / cfg.radar.microsat_antenna
off_nadir = np.radians(cfg.radar.off_nadir)
# calculate orbit
alt0 = cfg.orbit.alt
inc = cfg.orbit.inc
e = cfg.orbit.e
re_cycle_num = cfg.orbit.re_cycle_num 
re_day = cfg.orbit.re_day
if np.size(alt0)>1:
    alt = np.zeros_like(alt0)
    v = np.zeros_like(alt0)
    dis = np.zeros_like(alt0)
    for i in range(np.size(alt0)):
        alt[i], v[i], dis[i] = orbit_com(alt0[i], inc, e, 
                                           re_cycle_num, re_day, 
                                           cfg.orbit.con_repeat)
else:
    alt, v, dis = orbit_com(alt0, inc, e, re_cycle_num, re_day)
alt = alt / np.cos(a_bw / 2 + off_nadir)
  


B_amb = l0 * alt / 2 * np.sin(np.radians(cfg.scene.slope)) / baseline.reshape((baseline.shape+(1,)))


if np.size(cfg.scene.slope)>1:
    # antenna & peak power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(np.size(cfg.scene.slope)):
        plt.plot(baseline, B_amb[:,i])
        leg = np.append(leg, ["Scene slope=%d deg"%(cfg.scene.slope[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Baseline length (m)', fontsize=16)
    plt.ylabel('Height ambiguity (m)', fontsize=16)
#    plt.title(concept)
    plt.legend(leg)
#    plt.ylim(0, 2)
    plt.grid() 
    
if np.size(alt)>1:
    # antenna & peak power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(np.size(alt)):
        plt.plot(baseline, B_amb[:,i])
        leg = np.append(leg, ["Orbit height=%d km"%(cfg.orbit.alt[i]/1000)])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Baseline length (m)', fontsize=16)
    plt.ylabel('Height ambiguity (m)', fontsize=16)
#    plt.title(concept)
    plt.legend(leg)
#    plt.ylim(0, 2)
    plt.grid() 