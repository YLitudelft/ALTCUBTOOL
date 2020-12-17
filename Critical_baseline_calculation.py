# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 12:11:32 2019

This file is for calculating the critical baselines

@author: lyh
"""
# Setting up the packages
import numpy as np
from matplotlib import pyplot as plt
from constants import constants as const
from orbit.orbit import orbit_com
import ocs_io as tpio

# Setting up the paths
cfg_file = r"D:\research\TU Delft\programs\Alticube_tool\Alticubes'tool\parameters.cfg"
cfg = tpio.ConfigFile(cfg_file)

# lambda [m]
l0 = const.c / cfg.radar.f0
# along tr. bw [rad]
a_bw = l0 / cfg.radar.microsat_antenna
#off_nadir = np.radians(cfg.radar.off_nadir)
off_nadir = cfg.radar.off_nadir
off_nadir = np.radians(off_nadir)
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
#alt = alt / np.cos(a_bw / 2 + off_nadir) 
# off-nadir angle here is the smallest incidence angle
D = (np.tan(off_nadir) * alt + a_bw / 2 * alt) * 2 
off_nadir_max = np.arctan((D / 2 + off_nadir / 2 * alt) / alt)
incidence_angle = np.degrees((off_nadir_max + off_nadir) / 2)
if np.size(cfg.radar.bandwidth)>1:
    # critical baseline determined by the smallest incidence angle
    B_cri = 2 * cfg.radar.bandwidth.reshape((cfg.radar.bandwidth.shape+(1,))) * l0 * alt * np.tan(off_nadir) / const.c 
    cross_res  = 1 / np.sqrt(10**(np.cos(np.radians(incidence_angle))**3 * cfg.loss.SNR/ 10)) * l0 * alt / (2 * np.pi * B_cri)
    B_roll_e = 1e-3 / (B_cri/2) * a_bw * alt  
    n_r = np.sqrt((100 * alt / np.cos(np.radians(incidence_angle)) * np.sin(np.radians(incidence_angle)) * 
              l0 / 2 / np.pi / (B_cri/2) * np.sqrt((1 - gamma**2)/ gamma**2 / 2 / look_num))**2 * 0.267 /
              (1/(plu_res/1000)-1/(up_bound/1000)) * (1/(scale_az/1000)-1/(up_bound/1000)))
elif np.size(cfg.loss.SNR)>1:
    B_cri = 2 * cfg.radar.bandwidth * l0 * alt * np.tan(off_nadir) / const.c / np.cos(off_nadir)
    cross_res  = 1 / np.sqrt(10**(np.cos(np.radians(incidence_angle))**3 * cfg.loss.SNR.reshape((cfg.loss.SNR.shape+(1,)))/ 10)) * l0 * alt / (2 * np.pi * B_cri)
    B_roll_e = 1e-3 / B_cri * np.tan(off_nadir) * alt 

#if np.size(cfg.radar.bandwidth)>1:
#    # Critial baseline  
#    # plot figures
#    plt.figure()
#    plt.plot(cfg.radar.bandwidth/1e6, B_cri)
#    plt.xticks(fontsize=16)
#    plt.yticks(fontsize=16)
#    plt.xlabel('Bandwidth (MHz)', fontsize=16)
#    plt.ylabel('Critical baseline (m)', fontsize=16)
#
##    plt.ylim(0, 2)
#    plt.grid() 
#    # accross-track resolution  
#    # plot figures
#    plt.figure()
#    plt.plot(cfg.radar.bandwidth/1e6, cross_res)
#    plt.xticks(fontsize=16)
#    plt.yticks(fontsize=16)
#    plt.xlabel('Bandwidth (MHz)', fontsize=16)
#    plt.ylabel('Cross-track resolution (m)', fontsize=16)
#
##    plt.ylim(0, 2)
#    plt.grid() 

#if np.size(cfg.radar.microsat_antenna)>1:
#    # accross-track resolution  
#    # plot figures
#    plt.figure()
#    plt.plot(cfg.radar.microsat_antenna, cross_res)
#    plt.xticks(fontsize=16)
#    plt.yticks(fontsize=16)
#    plt.xlabel('Cross-track antenna dimension of the Microsat (m)', fontsize=16)
#    plt.ylabel('Cross-track resolution (m)', fontsize=16)
#
##    plt.ylim(0, 2)
#    plt.grid() 
# plot results
if np.size(cfg.radar.microsat_antenna)>1:
    # antenna & peak power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(np.size(a_bw)):
        plt.plot(cfg.radar.bandwidth/1e6, B_cri[:,i])
        leg = np.append(leg, ["Microsat's antenna=%.2f m"%(cfg.radar.microsat_antenna[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Bandwidth (MHz)', fontsize=16)
    plt.ylabel('Max. critical baseline (m)', fontsize=16)
#    plt.title(concept)
    plt.legend(leg)
#    plt.ylim(0, 2)
    plt.grid() 
    
if np.size(off_nadir)>1:
    off_nadir = np.degrees(off_nadir)
    if np.size(cfg.radar.bandwidth)>1:
        # antenna & peak power   
        # plot figures
        leg = []
        plt.figure()
        if np.size(cfg.radar.bandwidth)>1:
            for i in range(np.size(cfg.radar.bandwidth)):
                plt.plot(incidence_angle, B_cri[i,:])
                leg = np.append(leg, ["Bandwidth=%d MHz"%(cfg.radar.bandwidth[i]/ 1e6)])
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.xlabel('Incidence angle (deg)', fontsize=16)
            plt.ylabel('Critical baseline (m)', fontsize=16)
            #    plt.title(concept)
            plt.legend(leg)
            #    plt.ylim(0, 2)
            plt.grid() 
        else:
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            plt.plot(incidence_angle, B_cri)
            ax1.tick_params(labelsize=20)
            ax1.tick_params(labelsize=20)
            ax1.set_xlabel('Incidence angle (deg)', fontsize=20)
            ax1.set_ylabel('Critical baseline (m)', fontsize=20)
            ax2 = ax1.twinx()  
            plt.plot(incidence_angle, n_r,'r')
            ax2.set_ylabel('Altimeter noise (cm)', fontsize=20)
            ax2.tick_params(labelsize=20)
            ax2.tick_params(labelsize=20)
            fig.legend(['Critical baseline','Altimeter noise'],prop={'size': 20}, 
            loc=1, bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)        
        leg = []
        plt.figure()
        for i in range(np.size(cfg.radar.bandwidth)):
            plt.plot(incidence_angle, cross_res[i,:])
            leg = np.append(leg, ["Bandwidth=%d MHz"%(cfg.radar.bandwidth[i]/1e6)])
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('Incidence angle (deg)', fontsize=16)
        plt.ylabel('Cross-track resolution (m)', fontsize=16)
    #    plt.title(concept)
        plt.legend(leg)
    #    plt.ylim(0, 2)
        plt.grid() 
        
        plt.figure()
        plt.plot(incidence_angle,D/1000)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('Incidence angle (deg)', fontsize=16)
        plt.ylabel('Distance between the microsat and cubsat (km)', fontsize=16)
        plt.grid()
        plt.figure()
        plt.plot(D / 1000, off_nadir)
        plt.plot(D / 1000, np.degrees(off_nadir_max))
        plt.fill_between(D / 1000, off_nadir, np.degrees(off_nadir_max), facecolor = "yellow")
        plt.plot(D / 1000, incidence_angle, '--')
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('Distance between the active cubesat and the nearest passive CubSat (km)', fontsize=16)
        plt.ylabel('Incidence angle (deg)', fontsize=16)
        plt.grid()
        plt.legend(['Near range', 'Far range', 'Average'])
        # basline roll angle
        leg = []
        plt.figure()
        for i in range(np.size(cfg.radar.bandwidth)):
            plt.plot(incidence_angle, B_roll_e[i,:]*100)
            leg = np.append(leg, ["Bandwidth=%d MHz"%(cfg.radar.bandwidth[i]/1e6)])
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('Incidence angle (deg)', fontsize=16)
        plt.ylabel('Baseline roll error (cm)', fontsize=16)
    #    plt.title(concept)
        plt.legend(leg)
    #    plt.ylim(0, 2)
        plt.grid() 
    elif np.size(cfg.loss.SNR)>1:
        # antenna & peak power   
        # plot figure        
        leg = []
        plt.figure()
        for i in range(np.size(cfg.radar.off_nadir)):
            plt.plot(cfg.loss.SNR, cross_res[:,i])
            leg = np.append(leg, ["Incidence angle=%.1f deg"%(incidence_angle[i])])
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('SNR (dB)', fontsize=16)
        plt.ylabel('Cross-track resolution (m)', fontsize=16)
    #    plt.title(concept)
        plt.legend(leg)
    #    plt.ylim(0, 2)
        plt.grid() 
        
