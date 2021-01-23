# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 13:55:33 2019


     This file calculates the parameters for swath altimeters


@author: lyh
"""

# Setting up the packages
import numpy as np
from matplotlib import pyplot as plt
from constants import constants as const
from orbit.orbit import orbit_com
import ocs_io as tpio

# paths of the documents
cfg_file = r"D:\research\TU Delft\programs\ALTCUBTOOL\cfg\parameters.cfg"
cfg = tpio.ConfigFile(cfg_file)

# lambda [m]
l0 = const.c / cfg.radar.f0
# along tr. bw [rad]
a_bw = l0 / cfg.radar.microsat_antenna
# off-nadir angle if swath altimeter is implemented
off_nadir = np.radians(cfg.radar.off_nadir)
# calculate orbit
alt0 = cfg.orbit.alt
inc = cfg.orbit.inc
e = cfg.orbit.e
R = const.r_earth
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

# calculate the incidence angles at the near-range and the far-range
if np.size(off_nadir)==1:
    angle1 = np.arcsin(np.sin(off_nadir + a_bw.reshape((a_bw.shape+(1,))) / 2) * (alt + R) / R) - (off_nadir + a_bw.reshape((a_bw.shape+(1,))) / 2)
    angle2 = np.arcsin(np.sin(off_nadir - a_bw.reshape((a_bw.shape+(1,))) / 2) * (alt + R) / R) - (off_nadir - a_bw.reshape((a_bw.shape+(1,))) / 2)
else:
    angle1 = np.arcsin(np.sin(off_nadir.reshape((off_nadir.shape+(1,))) + a_bw / 2) * (alt + R) / R) - (off_nadir.reshape((off_nadir.shape+(1,))) + a_bw / 2)
    angle2 = np.arcsin(np.sin(off_nadir.reshape((off_nadir.shape+(1,))) - a_bw / 2) * (alt + R) / R) - (off_nadir.reshape((off_nadir.shape+(1,))) - a_bw / 2)
# calculate the swath length
swath_length = (angle1 - angle2) * R

if (np.size(off_nadir)>1 and np.size(alt)==1):
    # antenna & peak power   
    # plot figures   
    leg = []
    plt.figure()
    plt.plot(cfg.radar.off_nadir, swath_length / 1000)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Off nadir angle (deg)', fontsize=16)
    plt.ylabel('Swath length (km)', fontsize=16)
#    plt.title(concept)
#    plt.ylim(0, 2)
    plt.grid() 
if np.size(alt)>1:
    # plot figures
    leg = []
    plt.figure()
    if np.size(off_nadir)>1:
        for i in range(swath_length.shape[1]):
            plt.plot(cfg.radar.off_nadir, swath_length[:,i] / 1000)
            leg = np.append(leg, ["Orbit height=%d km"%(cfg.orbit.alt[i]/1000)])
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('Off nadir angle (deg)', fontsize=16)
        plt.ylabel('Swath length (km)', fontsize=16)
        plt.title(concept)
        plt.legend(leg)
    #    plt.ylim(0, 30)
        plt.grid()
    else:
        for i in range(swath_length.shape[1]):
            plt.plot(cfg.radar.microsat_antenna, swath_length[:,i] / 1000)
            leg = np.append(leg, ["Orbit height=%d km"%(cfg.orbit.alt[i]/1000)])
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('Cross track antenna size of the Microsat', fontsize=16)
        plt.ylabel('Swath length (km)', fontsize=16)
        plt.title(concept)
        plt.legend(leg)
    #    plt.ylim(0, 30)
        plt.grid()
        

    
