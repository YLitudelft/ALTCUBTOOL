# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:21:45 2020

This file deals with the performance analysis of the altimeters under different 
system parameters and concepts
Three concepts are considered here
1) Comb constellation
2) Specular reflection: microsat + cubesat
3) Swath interferometric altimeter: based on the SWOT's along-track SSH error spectrum 

@author: lyh
"""
# Setting up the packages
import numpy as np
from matplotlib import pyplot as plt
import ocs_io as tpio
from performance.performance import performance_tool

# Setting up paths
cfg_file = r"D:\research\TU Delft\programs\ALTCUBTOOL\cfg\parameters.cfg"
cfg = tpio.ConfigFile(cfg_file)
concept = 'Specular Constellation'##"#cfg.sim.concept'Comb Constellation'#
(P_r1, P_a1, ac1, raw_data_rate1, proc_rate1) = performance_tool(cfg_file, concept, 
                                                                 SAR_mode_enable=True, 
                                                                 interleaved_mode_enable=False)
###############################################################################
# ploting figure parts
if (np.size(cfg.radar.pulse_coh_bur)>1):
    # power
    plt.figure()
    plt.plot(cfg.radar.pulse_coh_bur,P_r1)
    plt.plot(cfg.radar.pulse_coh_bur,P_a1)
    leg = ["Peak power","Average power"]
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Coherent pulse number', fontsize=16)
    plt.ylabel('Required power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 0.5)
    plt.grid()
    # accuracy
    plt.figure()
    plt.plot(cfg.radar.pulse_coh_bur,ac1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Coherent pulse number', fontsize=16)
    plt.ylabel('Accuracy (cm)', fontsize=16)
    plt.title(concept)
    plt.ylim(0, 20)
    plt.grid()
    
    # raw data rate
    plt.figure()
    plt.plot(cfg.radar.pulse_coh_bur,raw_data_rate1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Coherent pulse number', fontsize=16)
    plt.ylabel('Raw data rate (Mb/s)', fontsize=16)
    plt.title(concept)
    plt.ylim(0, 3)
    plt.grid()
    
    # processed data rate
    plt.figure()
    plt.plot(cfg.radar.pulse_coh_bur,proc_rate1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Coherent pulse number', fontsize=16)
    plt.ylabel('Processed data rate (kb/s)', fontsize=16)
    plt.title(concept)
    plt.ylim(0, 50)
    plt.grid()
    
if (np.size(cfg.radar.bur_rep_int)>1):
    # power
    plt.figure()
    plt.plot(cfg.radar.bur_rep_int,P_a1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Burst interval time (ms)', fontsize=16)
    plt.ylabel('Required average power (W)', fontsize=16)
    plt.title(concept)
    plt.ylim(0, 0.06)
    plt.grid()
    # accuracy
    plt.figure()
    plt.plot(cfg.radar.bur_rep_int,ac1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Burst interval time (ms)', fontsize=16)
    plt.ylabel('Accuracy (cm)', fontsize=16)
    plt.title(concept)
    plt.ylim(0, 5)
    plt.grid()
    
    # raw data rate
    plt.figure()
    plt.plot(cfg.radar.bur_rep_int,raw_data_rate1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Burst interval time (ms)', fontsize=16)
    plt.ylabel('Raw data rate (Mb/s)', fontsize=16)
    plt.title(concept)
    plt.ylim(0, 20)
    plt.grid()
    
if (np.size(cfg.radar.bandwidth)>1):
    # power
    plt.figure()
    plt.plot(cfg.radar.bandwidth/1e6,P_r1)
    plt.plot(cfg.radar.bandwidth/1e6,P_a1)
    leg = ["Peak power","Average power"]
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Bandwidth (MHz)', fontsize=16)
    plt.ylabel('Required power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 0.5)
    plt.grid()
    # accuracy
    plt.figure()
    plt.plot(cfg.radar.bandwidth/1e6,ac1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Bandwidth (MHz)', fontsize=16)
    plt.ylabel('Accuracy (cm)', fontsize=16)
    plt.title(concept)
    plt.ylim(0, 5)
    plt.grid()
if (np.size(cfg.orbit.alt)>1):
    # power
    plt.figure()
    plt.plot(cfg.orbit.alt/1e3,P_r1)
    plt.plot(cfg.orbit.alt/1e3,P_a1)
    leg = ["Peak power","Average power"]
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Orbit height (km)', fontsize=16)
    plt.ylabel('Required power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 2)
    plt.grid()
    
if (np.size(cfg.radar.microsat_antenna_along)>1):
    # power
    plt.figure()
    plt.plot(cfg.radar.microsat_antenna_along,P_r1)
    plt.plot(cfg.radar.microsat_antenna_along,P_a1)
    leg = ["Peak power","Average power"]
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Along-track antenna dimension of the Microsat (m)', fontsize=16)
    plt.ylabel('Required power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 10)
    plt.grid()
    # accuracy
    plt.figure()
    plt.plot(cfg.radar.microsat_antenna_along,ac1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Cross-track antenna dimension of the Microsat (m)', fontsize=16)
    plt.ylabel('Altimeter noise (cm)', fontsize=16)
    plt.title(concept)
    #plt.ylim(3, 10)
    plt.grid()
if (np.size(cfg.radar.microsat_antenna)>1):
    # power
    plt.figure()
    plt.plot(cfg.radar.microsat_antenna,P_r1)
    plt.plot(cfg.radar.microsat_antenna,P_a1)
    leg = ["Peak power","Average power"]
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Cross-track antenna dimension of the Microsat (m)', fontsize=16)
    plt.ylabel('Required power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 200)
    plt.grid()
    # swath
    l0 = const.c / cfg.radar.f0
    plt.figure()
    plt.plot(cfg.radar.microsat_antenna,l0 / cfg.radar.microsat_antenna *cfg.orbit.alt/1000)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Cross-track antenna dimension of the Microsat (m)', fontsize=16)
    plt.ylabel('Swath (km)', fontsize=16)
    plt.title(concept)
    plt.grid()
if (np.size(cfg.radar.off_nadir)>1):
    # power
    plt.figure()
    plt.plot(cfg.radar.off_nadir,P_r1)
    plt.plot(cfg.radar.off_nadir,P_a1)
    leg = ["Peak power","Average power"]
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Incidence angle near range (deg)', fontsize=16)
    plt.ylabel('Required power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 1000)
    plt.grid()



    

