# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 17:11:18 2020

This file is the main function to conduct performance analysis of the altimeters under different concepts
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


cfg_file = r"D:\research\TU Delft\programs\Alticube_tool\Alticubes'tool\parameters.cfg"
cfg = tpio.ConfigFile(cfg_file)
concept = 'Specular Constellation'#'Comb Constellation'##"Specular Constellation"#cfg.sim.concept
(P_r, P_a, ac, raw_data_rate, proc_rate) = performance_tool(cfg_file, concept, 
                                                            SAR_mode_enable=False, 
                                                            interleaved_mode_enable=False)

###############################################################################
# ploting figure parts
if (np.size(cfg.radar.az_ant)>1) or (np.size(cfg.radar.microsat_antenna))>1:
    # antenna & peak power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(P_r.shape[1]):
        plt.plot(cfg.loss.SNR, P_r[:,i])
        if (np.size(cfg.radar.microsat_antenna))>1:
            leg = np.append(leg, ["Microsat's antenna=%.2f m"%(cfg.radar.microsat_antenna[i])])
        else:
            leg = np.append(leg, ["Antenna=%.2f m"%(cfg.radar.az_ant[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required peak power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 10)
    plt.grid()
    
    # antenna & average power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(P_a.shape[1]):
        plt.plot(cfg.loss.SNR, P_a[:,i])
        if (np.size(cfg.radar.microsat_antenna))>1:
            leg = np.append(leg, ["Microsat's antenna=%.2f m"%(cfg.radar.microsat_antenna[i])])
        else:
            leg = np.append(leg, ["Antenna=%.2f m"%(cfg.radar.az_ant[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required average power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 2)
    plt.grid()
    
    
    # antenna & accuracy   
    # plot figures
    leg = []
    plt.figure()
    for i in range(ac.shape[1]):
        plt.plot(cfg.loss.SNR, ac[:,i])
        if (np.size(cfg.radar.microsat_antenna))>1:
            leg = []
        else:
            leg = np.append(leg, ["Antenna=%.2f m"%(cfg.radar.az_ant[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Accuracy (cm)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
#    plt.ylim(2, 5)
    plt.grid()
    

# coherent burst num & accuracy   
# plot figures
if np.size(cfg.radar.pulse_coh_bur)>1:
    leg = []
    plt.figure()
    for i in range(P_r.shape[1]):
        plt.plot(cfg.loss.SNR, P_r[:,i])
        leg = np.append(leg, ["Pulses coh. burst num=%d"%(cfg.radar.pulse_coh_bur[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required peak power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 2)
    plt.grid()
    
    # antenna & average power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(P_a.shape[1]):
        plt.plot(cfg.loss.SNR, P_a[:,i])
        leg = np.append(leg, ["Pulses coh. burst num=%d"%(cfg.radar.pulse_coh_bur[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required average power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 0.5)
    plt.grid()
    leg = []
    plt.figure()
    for i in range(ac.shape[1]):
        plt.plot(cfg.loss.SNR, ac[:,i])
        leg = np.append(leg, ["Pulses coh. burst num=%d"%(cfg.radar.pulse_coh_bur[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Accuracy (cm)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(2, 8)
    plt.grid()


# burst rep interval & accuracy   
# plot figures
if np.size(cfg.radar.bur_rep_int)>1:
    leg = []
    plt.figure()
    for i in range(P_r.shape[1]):
        plt.plot(cfg.loss.SNR, P_r[:,i])
        leg = np.append(leg, ["burst rep interval=%d ms"%(cfg.radar.bur_rep_int[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required peak power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 10)
    plt.grid()
    
    # antenna & average power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(P_a.shape[1]):
        plt.plot(cfg.loss.SNR, P_a[:,i])
        leg = np.append(leg, ["burst rep interval=%d ms"%(cfg.radar.bur_rep_int[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required average power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 0.5)
    plt.grid()
    leg = []
    plt.figure()
    for i in range(ac.shape[1]):
        plt.plot(cfg.loss.SNR, ac[:,i])
        leg = np.append(leg, ["burst rep interval=%d ms"%(cfg.radar.bur_rep_int[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Accuracy (cm)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(1, 6)
    plt.grid()
    
    
if np.size(cfg.radar.bandwidth)>1:    
    # bandwidth & accuracy  
    leg = []
    plt.figure()
    for i in range(ac.shape[1]):
        plt.plot(cfg.loss.SNR, P_r[:,i])
        leg = np.append(leg, ["Bandwidth=%d MHz"%(cfg.radar.bandwidth[i]/1e6)])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required peak power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 1)
    plt.grid()
    
    # antenna & average power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(ac.shape[1]):
        plt.plot(cfg.loss.SNR, P_a[:,i])
        leg = np.append(leg, ["Bandwidth=%d MHz"%(cfg.radar.bandwidth[i]/1e6)])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required average power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 0.1)
    plt.grid()
    # plot figures
    leg = []
    plt.figure()
    for i in range(ac.shape[1]):
        plt.plot(cfg.loss.SNR, ac[:,i])
        leg = np.append(leg, ["Bandwidth=%d MHz"%(cfg.radar.bandwidth[i]/1e6)])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Accuracy (cm)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.grid()
    
if np.size(cfg.orbit.alt)>1:    
    # orbit & accuracy   
    # plot figures
    # antenna & peak power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(P_r.shape[1]):
        plt.plot(cfg.loss.SNR, P_r[:,i])
        leg = np.append(leg, ["Orbit height=%d km"%(cfg.orbit.alt[i]/1000)])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required peak power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 2)
    plt.grid()
    
    # antenna & average power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(P_a.shape[1]):
        plt.plot(cfg.loss.SNR, P_a[:,i])
        leg = np.append(leg, ["Orbit height=%d km"%(cfg.orbit.alt[i]/1000)])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required average power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 0.5)
    plt.grid()
    leg = []
    plt.figure()
    for i in range(ac.shape[1]):
        plt.plot(cfg.loss.SNR, ac[:,i])
        leg = np.append(leg, ["Orbit height=%d km"%(cfg.orbit.alt[i]/1000)])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Accuracy (cm)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(2, 8)
    plt.grid()
    
    
    # SNR & interferometric accuracy   
    # plot figures
    if cfg.sim.interferometer_enable:
        leg = []
        plt.figure()
        for i in range(ac.shape[1]):
            plt.plot(cfg.loss.SNR, cross_res[:,i])
            leg = np.append(leg, ["Orbit height=%d km"%(cfg.orbit.alt[i]/1000)])
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('SNR (dB)', fontsize=16)
        plt.ylabel('Cross-track resolution (m)', fontsize=16)
        plt.title(concept)
        plt.legend(leg)
    #    plt.ylim(0, 0.5)
        plt.grid()
#    leg = []
#    plt.figure()
#    for i in range(ac.shape[1]):
#        plt.plot(cfg.loss.SNR, height_ac[:,i] * 100)
#        leg = np.append(leg, ["Orbit height=%d km"%(cfg.orbit.alt[i]/1000)])
#    plt.xticks(fontsize=16)
#    plt.yticks(fontsize=16)
#    plt.xlabel('SNR (dB)', fontsize=16)
#    plt.ylabel('Height accuracy (cm)', fontsize=16)
#    plt.title(concept)
#    plt.legend(leg)
##    plt.ylim(0, 0.5)
#    plt.grid()
    
if (np.size(cfg.radar.off_nadir)>1):
    # antenna & peak power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(P_r.shape[1]):
        plt.plot(cfg.loss.SNR, P_r[:,i])
        leg = np.append(leg, ["off nadir angle=%.2f deg"%(cfg.radar.off_nadir[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required peak power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 200)
    plt.grid()
    
    # antenna & average power   
    # plot figures
    leg = []
    plt.figure()
    for i in range(P_a.shape[1]):
        plt.plot(cfg.loss.SNR, P_a[:,i])
        leg = np.append(leg, ["off nadir angle=%.2f deg"%(cfg.radar.off_nadir[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Required average power (W)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.ylim(0, 100)
    plt.grid()
if np.size(cfg.radar.baseline)>1:    
    # SNR & interferometric accuracy   
    # plot figures
    leg = []
    plt.figure()
    for i in range(cfg.radar.baseline.shape[0]):
        plt.plot(cfg.loss.SNR, cross_res[:,i])
        leg = np.append(leg, ["Baseline length=%.1f m"%(cfg.radar.baseline[i])])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('SNR (dB)', fontsize=16)
    plt.ylabel('Cross-track resolution (m)', fontsize=16)
    plt.title(concept)
    plt.legend(leg)
    plt.grid()