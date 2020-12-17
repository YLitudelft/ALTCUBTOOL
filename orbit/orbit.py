# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 11:13:05 2019

@author: lyh

    **Calculating orbit parameters**

     :alt0:        Initial orbit height [km]
     :inc:         Inclination [deg]
     :e:           Eccentricity
     :re_cycle_num:         Repeat cycle number [number/cycle]
     :alt:       Final orbit height [km]
     :v:         Satellite velocity [km/s]
     :dis:          Sampling distance at equator [km]

 """

import numpy as np
from constants import constants as const


def orbit_com(alt0, inc, e, re_cycle_num, re_day, con_repeat=False):
    
    threshold = 10
    max_in = 100
    i = 0
    if con_repeat:
        while i<max_in:
            G = const.g_star
            M = const.m_earth
            Mu = G * M
            orb_r = alt0 + const.r_earth
            sv = np.sqrt(Mu / orb_r)
            period = 2 * np.pi * orb_r / sv
            # Precession [deg/min]
            p = np.degrees(-3 / 2 * const.j2 * np.sqrt(Mu / orb_r**3) * (const.r_earth / orb_r)**2 *
                           np.cos(np.radians(inc)) / (1 - e**2)**2) * 60
            # Earth rotation [deg/min]
            earth_r = 360 / const.sd
            # Peprime
            Pe = 360 / (earth_r - p)
            
            P = re_day * Pe / re_cycle_num 
            v = (2 * np.pi * G * M / (P * 60))**(1/3)
            alt = v * 60 * P / 2 / np.pi - const.r_earth
            dis = 2 * np.pi * const.r_earth  / re_cycle_num
            if np.abs(alt - alt0)<threshold:
                i=100
            else:
                i=i+1
                alt0 = alt
    else:
        G = const.g_star
        M = const.m_earth
        Mu = G * M
        orb_r = alt0 + const.r_earth
        v = np.sqrt(Mu / orb_r)
        alt = alt0
        dis = 0   
    return alt, v, dis


