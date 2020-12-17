# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 13:40:23 2019

It is the key functions for performance analysis

@author: lyh
"""
import numpy as np
from matplotlib import pyplot as plt
from constants import constants as const
import ocs_io as tpio
from orbit.orbit import orbit_com
from performance.cal_functions import cell_calculate, area_cal

# THIS FUNCTION OFFERS THE PERFORMANCE ANALYSIS FOR CUBESAT SYSTEMS        
def performance_tool(cfg_file, concept, SAR_mode_enable, interleaved_mode_enable):
    # INITIAL PARAMETERS
    cfg = tpio.ConfigFile(cfg_file)
    # lambda [m]
    l0 = const.c / cfg.radar.f0
    # along tr. bw [rad]
    a_bw = l0 / cfg.radar.az_ant
    # cross tr. bw [rad]
    c_bw = l0 / cfg.radar.c_ant
    # along tr. eff. or smallest bw [deg]
    eff_a_bw = a_bw / np.sqrt(2)
    # cross tr. eff. or smallest bw [deg]
    eff_c_bw = c_bw / np.sqrt(2)
    # antenna gain
    gain_t = 10 * np.log10(4 * np.pi / a_bw / c_bw / 2)
    gain_r = gain_t
    # compressed pulse length
    pw = const.c / 2 / cfg.radar.bandwidth
    
    
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
    
    # mode recognize
    if concept=="Specular Constellation":
        # along tr. bw [rad]
        a_bw = l0 / cfg.radar.microsat_antenna_along
        # cross tr. bw [rad]
        c_bw = l0 / cfg.radar.microsat_antenna
        eff_a_bw = a_bw / np.sqrt(2)
        alt_m = alt / np.cos(a_bw / 2)
        gain_t = 10 * np.log10(4 * np.pi / a_bw / c_bw / 2) - 3
        if cfg.sim.max_coverage:
            gain_r = 10 * np.log10(4 * np.pi / c_bw / c_bw / 2)
        if cfg.sim.interferometer_enable:
            # based on Ulaby model if consider the ocean surface is mirror alph=1
            # at far range the SNR should be enough
            D = (np.tan(np.radians(cfg.radar.off_nadir)) * alt + c_bw / 2 * alt) * 2 
            off_nadir = np.arctan((D / 2 + np.radians(cfg.radar.off_nadir) / 2 * alt) / alt)# maximal off-nadir angle
            gain_r = gain_r - 10 * np.log10(4.34 * np.degrees(off_nadir)) - 3
        else:
             gain_r = gain_r - 3
    if SAR_mode_enable:
        pulse_coh_bur = cfg.radar.pulse_coh_bur
    else:
        pulse_coh_bur = 1
        
    
    if cfg.sim.plot_spec_coverage:
        if np.size(cfg.radar.microsat_antenna)>1:
            ang = l0 / cfg.radar.microsat_antenna 
            plt.figure()
            for i in range(np.size(alt_m)):
                plt.plot(cfg.radar.microsat_antenna, ang * cfg.orbit.alt[i]/1000)
                leg = np.append(leg, ["Orbit height=%d km"%(cfg.orbit.alt[i]/1000)])
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.xlabel('Cross-track size of microcube (m)', fontsize=16)
            plt.ylabel('Coverage area (km)', fontsize=16)
            plt.title(concept)
            plt.legend(leg)
        #        plt.ylim(2, 5)
            plt.grid()
    
    
    # calculation performance
    # Doppler BW [Hz] and prf
    Dop_bw = 2 * v / (l0 / eff_a_bw)
    if SAR_mode_enable:
        prf = 18700
    else:
        prf = 4000#int(1.4 * Dop_bw)
    prf_r = 1.4 * Dop_bw
    if cfg.sim.plot_prf:
        #plot the prf rule
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(l0 / a_bw, prf_r)
        #ax1.set_xticks(fontsize=16)
        #ax1.set_yticks(fontsize=16)
        ax1.set_xlabel('Along-track antenna dimension of the Microsat (m)', fontsize=16)
        ax1.set_ylabel('PRF (Hz)', fontsize=16, color="blue")
        ax1.set_title("PRF vs. swath")
        ax2 = ax1.twinx()  
        ax2.plot(l0 / a_bw, l0 / (1 / cfg.radar.microsat_antenna_along) * alt_m / 1000, 'r')
        ax2.set_ylabel('Cross-track swath length (km)', fontsize=16, color="red")
        plt.show()
        plt.grid()
    
    
    prt = 1 / prf
    # free space between pulses [m]
    #F_length = (prt - cfg.radar.pw_uc / 1000) / 2 * const.c
    # max # pulses in burst [ ]
    #max_pulse_burst = (alt * 2 / const.c / prt).astype(int)
    # along tr. footprint [km]
    al_res = eff_a_bw * alt
    # cross tr. footprint [km]
    # pulse interval time
    if SAR_mode_enable:
        bur_rep_int = cfg.radar.bur_rep_int
    else:
        bur_rep_int = prt * 1000
    c_res = eff_c_bw * alt
    # distance between bursts [m]
    dis_burst = bur_rep_int * v / 1000
    
    
    if cfg.sim.pulse_limited_enable:
        # pulse limited radius on ground [m]
        pul_rad = np.sqrt(2 * alt * pw)
        # Pulse lim. resolution [m]
        plu_res = 2 * pul_rad
    
    if SAR_mode_enable:
        # Synthetic Aperture length [m]
        syn_length = 1 / prf * pulse_coh_bur * v
        # SAR resolution [m]
        sar_res = prf * alt * l0 / (2 * pulse_coh_bur * v)
        # Raney limit [m]
        raney_l = 3 * pw * alt / sar_res
        # cells in Raney region
        cel = cell_calculate(al_res, raney_l, sar_res)
        # burst length [s]
        bur_length = pulse_coh_bur * prt
        # incoh. observ. per second
        N_incoh = cel / bur_rep_int * 1000
        # incoh. observ. per res. cell
        N_incoh_r = cel * sar_res / v / bur_rep_int * 1000
    else:
        # cells in Raney region
        cel = 1
        sar_res = plu_res
        # burst length [s]
        bur_length = cfg.radar.pw_uc
        if np.size(dis_burst)>1 and np.size(sar_res)>1 :
            N_incoh = np.zeros_like(dis_burst)
            for i in range(np.size(N_incoh)):
                if dis_burst[i]>0.305 * alt[i] * l0 / sar_res[i]:
                    N_incoh[i] = cel / bur_rep_int * 1000
                else:
                    N_incoh[i] = v / (0.305 * alt[i] * l0 / sar_res[i])
        elif np.size(dis_burst)>1 and np.size(sar_res)==1 :
            N_incoh = np.zeros_like(dis_burst)
            for i in range(np.size(N_incoh)):
                if dis_burst[i]>0.305 * alt * l0 / sar_res:
                    N_incoh[i] = cel / bur_rep_int[i] * 1000
                else:
                    N_incoh[i] = v / (0.305 * alt * l0 / sar_res)
        elif np.size(sar_res)>1:
            N_incoh = np.zeros_like(sar_res)
            for i in range(np.size(sar_res)):
                if dis_burst>0.305 * alt * l0 / sar_res[i]:
                    N_incoh[i] = cel / bur_rep_int * 1000
                else:
                    N_incoh[i] = v / (0.305 * alt[i] * l0 / sar_res)
        else:
            if dis_burst>0.305 * alt * l0 / sar_res:
                N_incoh = cel / bur_rep_int * 1000
            else:
                N_incoh = v / (0.305 * alt * l0 / sar_res)
        N_incoh_r = N_incoh * sar_res / v
        
        
    # power performance
    KTF = -204 + cfg.loss.NF
    # illuminated area [m2]
    A = area_cal(sar_res, plu_res, pul_rad)
    A0 = area_cal(plu_res, plu_res, pul_rad)
    
    
    # Power calculation & accuracy
    sigma = 10 * np.log10(A) + cfg.loss.sigma0
    sigma_w = 10 * np.log10(A0) + cfg.loss.sigma0
    # req. peak power [W]
    if (np.size(cfg.loss.SNR)>1):
        P_r = ((4 * np.pi)**3 * alt**4 * 10**((cfg.loss.SNR.reshape((cfg.loss.SNR.shape+(1,))) + cfg.loss.loss_in_at + KTF - gain_t - gain_r - sigma) / 10) / 
                l0**2 / (cfg.radar.pw_uc / 1000) / pulse_coh_bur)
        # req. av. power [W]
        P_a = P_r * pulse_coh_bur * cfg.radar.pw_uc / bur_rep_int
    else:
        P_r = ((4 * np.pi)**3 * alt**4 * 10**((cfg.loss.SNR + cfg.loss.loss_in_at + KTF - gain_t - gain_r - sigma) / 10) / 
            l0**2 / (cfg.radar.pw_uc / 1000) / pulse_coh_bur)
        # req. av. power [W]
        P_a = P_r * pulse_coh_bur * cfg.radar.pw_uc / bur_rep_int
    
    # calculating the SNR before the compression processing
    P_rc = cfg.radar.input_power * 10**((gain_t + gain_r - cfg.loss.loss_in_at + sigma_w) / 10) * l0**2 / (4 * np.pi)**3 / alt**4 
    noise = 10**(KTF/10) * cfg.radar.bandwidth
    SNR0 = 10 * np.log10(P_rc / noise)
    P_rcw = 10 * np.log10(P_rc) + 30
    noise_w = 10 * np.log10(noise) + 30
    
    # single cell height accuracy [cm]
    if (np.size(cfg.loss.SNR)>1):
        ac = pw / np.sqrt(2 * N_incoh_r * (1 / (1 + 1/10**(cfg.loss.SNR.reshape((cfg.loss.SNR.shape+(1,))) / 10)))) * 100
    else:
        ac = pw / np.sqrt(2 * N_incoh_r * (1 / (1 + 1/10**(cfg.loss.SNR / 10)))) * 100
    
    # raw data rate [Mb/s], 16 bit samples
    raw_data_rate = pulse_coh_bur * cfg.radar.rg_num * 1000 / bur_rep_int / 1e6 * 16
    # proc.per channel [kb/s]
    if SAR_mode_enable:
        proc_rate = v / sar_res * cfg.radar.rg_num * 16 / 1000
    else:
        proc_rate = 4 * cfg.radar.rg_num * 16 / 1000
    
    if cfg.sim.interferometer_enable:
        # cross track resolution [m]
        if np.size(cfg.loss.SNR)>1:
            cross_res  = 1 / np.sqrt(10**(np.cos(a_bw / 2)**3 * cfg.loss.SNR.reshape((cfg.loss.SNR.shape+(1,))) / 10)) * l0 * alt / (2 * np.pi * cfg.radar.baseline)
    # This calculation is meaningless
    #    height_ac = cross_res / np.cos(a_bw / 2)
        
        
    if interleaved_mode_enable:
        P_a = P_r * prf *  cfg.radar.pw_uc / 1000
        if SAR_mode_enable:
            N_incoh = cel / (pulse_coh_bur / prf)
        else:
        # cells in Raney region
            if 1/prf * v>0.305 * alt * l0 / sar_res:
                N_incoh = cel / pulse_coh_bur / prf
            else:
                N_incoh = min(v / (0.305 * alt * l0 / sar_res),prf)
        # incoh. observ. per res. cell
        N_incoh_r = cel * sar_res / v / pulse_coh_bur *prf
        ac = pw / np.sqrt(2 * N_incoh_r * (1 / (1 + 1/10**(cfg.loss.SNR.reshape((cfg.loss.SNR.shape+(1,))) / 10)))) * 100
        raw_data_rate = prf * cfg.radar.rg_num / 1e6 * 16
        
        # analysis of cells in the raney region
    if cfg.sim.plot_raney:
        plt.figure()
        plt.plot(l0 / a_bw, 2*raney_l/1000)
        plt.xlabel('Along-track antenna dimension of the Microsat (m)', fontsize=16)
        plt.ylabel('Along-track footprint  (km)', fontsize=16)
        plt.plot(l0 / a_bw, al_res/1000, 'r')
        ax2.set_ylabel('Length (km)', fontsize=16)
        plt.legend(['Raney region length','Along-track footprint'])
        plt.show()
        plt.grid()
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(l0 / a_bw, cel)
        #ax1.set_xticks(fontsize=16)
        #ax1.set_yticks(fontsize=16)
        ax1.set_xlabel('Along-track antenna dimension of the Microsat (m)', fontsize=16)
        ax1.set_ylabel('cell number in the raney region', fontsize=16, color="blue")
        ax2 = ax1.twinx()  
        ax2.plot(l0 / a_bw, ac, 'r')
        ax2.set_ylabel('SSH accuracy', fontsize=16, color="red")
        plt.show()
        plt.grid()
        
    return P_r, P_a, ac, raw_data_rate, proc_rate


    




    
    



