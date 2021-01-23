# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 17:19:31 2020

This file focuses on the performance under different modes

@author: lyh
"""
# Setting up the packages
import numpy as np
from matplotlib import pyplot as plt
import ocs_io as tpio
from performance import performance_tool

# Setting up paths
cfg_file = r"D:\research\TU Delft\programs\ALTCUBTOOL\cfg\parameters.cfg"
cfg = tpio.ConfigFile(cfg_file)


concept = 'Comb Constellation'#"Specular Constellation"#cfg.sim.concept
(P_r1, P_a1, ac1, raw_data_rate1, proc_rate1) = performance_tool(cfg_file, concept, 
                                                                 SAR_mode_enable=False, 
                                                                 interleaved_mode_enable=False)
(P_r2, P_a2, ac2, raw_data_rate2, proc_rate2) = performance_tool(cfg_file, concept, 
                                                                 SAR_mode_enable=True, 
                                                                 interleaved_mode_enable=False)
(P_r3, P_a3, ac3, raw_data_rate3, proc_rate3) = performance_tool(cfg_file, concept, 
                                                                 SAR_mode_enable=True, 
                                                                 interleaved_mode_enable=True)
###############################################################################
# ploting figure parts
# peak power
plt.figure()
plt.plot(cfg.loss.SNR,P_r1)
plt.plot(cfg.loss.SNR,P_r2)
plt.plot(cfg.loss.SNR,P_r3)
leg = ["Traditional mode","Burst SAR mode","Interleaved SAR mode"]
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('SNR (dB)', fontsize=16)
plt.ylabel('Required peak power (W)', fontsize=16)
plt.title(concept)
plt.legend(leg)
plt.ylim(0, 5)
plt.grid()

#average power
plt.figure()
plt.plot(cfg.loss.SNR,P_a1)
plt.plot(cfg.loss.SNR,P_a2)
plt.plot(cfg.loss.SNR,P_a3)
leg = ["Traditional mode","Burst SAR mode","Interleaved SAR mode"]
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('SNR (dB)', fontsize=16)
plt.ylabel('Required average power (W)', fontsize=16)
plt.title(concept)
plt.legend(leg)
plt.ylim(0, 0.8)
plt.grid()

# accuracy
plt.figure()
plt.plot(cfg.loss.SNR,ac1)
plt.plot(cfg.loss.SNR,ac2)
plt.plot(cfg.loss.SNR,ac3)
leg = ["Traditional mode","Burst SAR mode","Interleaved SAR mode"]
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('SNR (dB)', fontsize=16)
plt.ylabel('Accuracy (cm)', fontsize=16)
plt.title(concept)
plt.legend(leg)
plt.ylim(0, 6)
plt.grid()

    