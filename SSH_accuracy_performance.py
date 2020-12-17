# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:31:57 2020

@author: lyh
It is the program to calculate the absolute SSH and relative SSH accuracy 
Three concepts are considered here
1) Comb constellation
2) Specular reflection: microsat + cubesat
3) Swath interferometric altimeter: based on the SWOT's along-track SSH error spectrum 
"""

# Setting up the packages
import numpy as np
from matplotlib import pyplot as plt
from constants import constants as const
import ocs_io as tpio
from error_functions import re_acuracy, re_acuracy_data,re_wet_one_radiometer, timing_error_data
from error_functions import re_acuracy_orbit_data,altimeter_rand, altimeter_rand_insar, re_accuracy_orbit_slope

# Setting up paths
cfg_file = r"D:\research\TU Delft\programs\Alticube_tool\Alticubes'tool\parameters.cfg"
cfg = tpio.ConfigFile(cfg_file)

# Input parameters
l0 = const.c / cfg.radar.f0
# compressed pulse length
pw = const.c / 2 / cfg.radar.bandwidth
# pulse limited radius on ground [m]
pul_rad = np.sqrt(2 * cfg.orbit.alt * pw)
# Pulse lim. resolution [m]
plu_res = 2 * pul_rad
# Weather parameters. It is used for calculate absolute SSB error
SWH = 2 #[m] 
# Antenna parameters
ant = 0.5 #[m]
#dis = l0 / ant * cfg.orbit.alt
D = 35000 # [km] distance between the microsat and the cubesat for interferometric case
B = 100 # [m] interferometric basline
Bpe = 1e-3 # [m] baseline roll accuracy
SNR = 11 # [dB]
gamma = 0.3 # coherence of the interfeograms from SWOT

# For the cubesat and microsat, the input parameters are
# antenna length
D_micro_c = 0.5 #[m]
D_micro_al = 1  #[m]
foot_print  = l0 / D_micro_c * cfg.orbit.alt # calculating resolution
# CALCULATING SYSTEM GEOMETRY
a_bw = l0 / ant
off_nadir_near = np.arctan((D / 2 - a_bw / 2 * cfg.orbit.alt) / cfg.orbit.alt)
off_nadir_max = np.arctan((D / 2 + a_bw / 2 * cfg.orbit.alt) / cfg.orbit.alt)
inc = np.degrees((off_nadir_max + off_nadir_near) / 2) # average off_nadir angle
res_insar = pw / np.sin(np.radians(inc)) # range resolution in interferometric case

# azimuth resolving scale for alticubes interferometric concept
scale_az = 1000#np.linspace(1000, 50000, 20)#5000
up_bound = 1000000 # 1000 km, a fixed value
# azimuth resolving scale for alticubes concept1 and 2 
# It is always evaluated by 1s averaging in along-track direction (refer to altika/crosat-2 manuscript)
scale_az1 = 7500#np.linspace(1000, 50000, 20)#5000

## pixel spacing in SWOT system
#pixel_SWOT = 7500 # cross-track averaging to 7.5km in SWOT
#Swath_swot = 50000 # 10-60 km

# Factor that link alticube and SWOT for the cross-track direction
# The errror spectrum given in SWOT has already been converted to white noise spectrum 
along_pixel = plu_res
# The gain of the system compared with the traitional nadir-looking altimeters
factor =  np.sqrt(along_pixel / plu_res) 
look_num = (1000 / 0.5) * (1000 / res_insar)
factor3 = np.sqrt(22250 / look_num) # 22250 from the SWOT possible averaging swath looks from 3700 to 40800
#factor_swot = 10 / 2 / 891000 / np.tan(np.radians(2.35)) * np.sqrt(22250)
# averging in along-track
N = scale_az / plu_res

# Accuracy initialization, it is mainly for the absolute accuracy
dis = 6340#50000
# The following power spectral indexes are from SWOT
p_ssb = -8/3
p_dr = -3
p_dw = -8/3
p1 = -1.79
p2 = -0.814
pi= -2.1

# axis for plot
Nx = 401
Ny = 301
gridx = np.linspace(-2200000, 2200000, Nx)
gridy = np.linspace(-1650000, 1650000, Ny)
#plot sign
ssbplot = False
ddplot = False
dwplot = False
diplot = False
alti_noise_plot = False
orbitplot = False


# THIS IS THE PART BEGAIN TO CALCULATE THE ERROR COMPONENTS
# Sources 1: sea state bias
# relative ssh for concept 1 & 2 from 10 km distance between small satellites 
# based on -3.8 power spectrum provided by SWOT
SSB_ab = 0.01 * SWH *100 / factor# 1% of the SWH at Ka band
SSB_r1 = np.sqrt(re_acuracy(dis/1000, SSB_ab, p_ssb))#0.78 / factor
#SSB_r2 = np.sqrt(re_acuracy(dis/2/1000,SSB_ab, p_ssb))#0.5 / factor
# for concept 3
SSB_r = np.sqrt(factor3 / 3 * 1e-9 / 2.8 * 
              ((up_bound/1000)**2.8 - (scale_az / 1000)**2.8))
SSB_data_comb = re_acuracy_data(SSB_ab*np.sqrt(2),p_ssb, Nx, Ny) # np.sqrt(2) is used to generate the relative error for concept 124
SSB_data_insar = re_acuracy_data(SSB_r,p_ssb, Nx, Ny)

# Source 2: dry tropospheric delay
dd_ab = 0.7 / factor
atp_a = 0.7 / factor #It is a semi-variance of mbar. A correction based on some external measurements
dd_abr = np.sqrt(atp_a * 2) * 2.277 /10
dd_r = np.sqrt(factor3 * 5 * 1e-9 / 2 * ((up_bound/1000)**2 - (scale_az / 1000)**2))
dd_r1 = dd_abr
dd_rswot=np.sqrt(re_acuracy(dis/1000,dd_ab, p_dr))
dd_data_comb = re_acuracy_data(dd_ab*np.sqrt(2),p_dr, Nx, Ny)
dd_data_insar = re_acuracy_data(dd_r,p_dr, Nx, Ny)

# Source 3: wet delay
# They are the functions from SWOT. It is the correction after 1-beam radiometer.
dw_ab = 1.2
dw_r1 = 0.4
if (np.size(scale_az)==1):
    if scale_az>1000/0.0023:
        dw_r = np.sqrt((9.5e-5 / 0.79 * ((up_bound/1000)**0.79 - (scale_az / 1000)**0.79)) * factor3)
    elif scale_az>1000/0.0683: 
        dw_r = np.sqrt((9.5e-5 / 0.79 * ((up_bound/1000)**0.79 - (1/0.0023)**0.79)+
                       0.036 / 0.186 * ((1/(scale_az / 1000))**0.186 - (0.0023)**0.186)) * factor3)
    else:
     dw_r = np.sqrt((9.5e-5 / 0.79 * ((up_bound/1000)**0.79 - (1/0.0023)**0.79)+
                       0.036 / 0.186 * ((0.0683)**0.186 - (0.0023)**0.186)+
                       0.32*(1 / (scale_az/1000) - 0.0683)) * factor3)
else:
    for i in range(np.size(scale_az)):
        if scale_az[i]>1000/0.0023:
            dw_r = np.sqrt((9.5e-5 / 0.79 * ((up_bound/1000)**0.79 - (scale_az[i] / 1000)**0.79)) * factor3)
        elif scale_az[i]>1000/0.0683: 
            dw_r = np.sqrt((9.5e-5 / 0.79 * ((up_bound/1000)**0.79 - (1/0.0023)**0.79)+
                       0.036 / 0.186 * ((1/(scale_az[i] / 1000))**0.186 - (0.0023)**0.186)) * factor3)
        else:
            dw_r = np.sqrt((9.5e-5 / 0.79 * ((up_bound/1000)**0.79 - (1/0.0023)**0.79)+
                       0.036 / 0.186 * ((0.0683)**0.186 - (0.0023)**0.186)+
                       0.32*(1 / (scale_az[i]/1000) - 0.0683)) * factor3)
dw_rswot=np.sqrt(re_wet_one_radiometer(dis/1000, dw_ab, p1, p2))
dw_data_comb = re_acuracy_data(dw_ab*np.sqrt(2),p_dw, Nx, Ny)
dw_data_insar = re_acuracy_data(dw_r,p_dr, Nx, Ny)
   
# Source 4: ionospheric delay
di_ab = 0.3 # [m] After the correction based on the Ionex model
di_abr = 2 * 40.28 / cfg.radar.f0**2 * 0.1e16 * 100 # It is corrected based on TEC measurements
di_r = np.sqrt(factor3 * 1e-8 / 1.1 * ((up_bound/1000)**1.1 - (scale_az / 1000)**1.1))
di_r1 = di_abr
di_rswot=np.sqrt(re_acuracy(dis/1000, di_ab, pi))
di_data_comb = re_acuracy_data(di_ab*np.sqrt(2),pi, Nx, Ny)
di_data_insar = re_acuracy_data(di_r,pi, Nx, Ny)

# Source 5: altimeter noise
n_ab = 1 / (scale_az1 / plu_res) #pw / np.sqrt(2 * (1 / (1 + 1 / 10**(SNR / 10)))) * 100
n_r1 = n_ab * np.sqrt(2)
n_r2 = n_ab * np.sqrt(2)
n_r = np.sqrt((100 * cfg.orbit.alt / np.cos(np.radians(inc)) * np.sin(np.radians(inc)) * 
              l0 / 2 / np.pi / B * np.sqrt((1 - gamma**2)/ gamma**2 / 2 / look_num))**2 * 0.267 /
              (1/(plu_res/1000)-1/(up_bound/1000)) * (1/(scale_az/1000)-1/(up_bound/1000)))

n_data_comb = altimeter_rand(Nx, Ny, n_r1)
n_data_insar = altimeter_rand_insar(Nx, Ny)
n_data_insar = n_data_insar / np.std(n_data_insar) * n_r

# Source 6: orbit error
n_orb_ab = 2 # These values have already been provided by the 1s averaging in SWOT/Altika documents
n_orb_abr = 0.09
# basline roll error for concept 3
kf = Bpe / B * cfg.orbit.alt * 100
n_orb_r = (off_nadir_max - off_nadir_near) * kf / 2  #np.sqrt(factor3 * 1.9631 / 0.9922 * 1e-5 * ((up_bound/1000)**0.9922 - (scale_az / 1000)**0.9922)) 
n_orb_r1 = n_orb_abr
n_orb_r2 = n_orb_abr * np.sqrt(2) / 2
n_orb_data_comb = re_acuracy_orbit_data(n_orb_r1, gridx, gridy)
n_orb_data_specular = re_acuracy_orbit_data(n_orb_r2, gridx, gridy)
grid_swath = np.linspace(off_nadir_near, off_nadir_max, 101)
n_orb_data_insar = re_accuracy_orbit_slope(gridx, gridy, kf)

# Source 7: timing error
st_os = 1e-10 # Depending on the stability of the oscillator  
st_syn = 1e-10  # Depending on the stability of the oscillator 
n_syn = (st_syn / 2) * const.c * 100 / (scale_az1 / plu_res)
n_t_ab1 = (st_os / 2) * const.c * 100 / (scale_az1 / plu_res)
n_t_ab2 = np.sqrt(n_t_ab1**2 + n_syn**2)
n_t_r1 = n_t_ab1 * np.sqrt(2)
n_t_r21 = n_t_ab1 / np.sqrt(2)
n_t_r2 = np.sqrt(n_t_r21**2 + n_syn**2)
n_t_r = np.sqrt(n_t_ab2**2 / look_num / (1/(plu_res/1000)-1/(up_bound/1000)) * (1/(scale_az/1000)-1/(up_bound/1000)))
n_timing_error_comb = timing_error_data(Nx, Ny, n_t_r1)
n_timing_error_specular = timing_error_data(Nx, Ny, n_t_r2)
n_timing_error_insar = timing_error_data(Nx, Ny, n_t_r)

# total absolute SSH error
T_ssh_ab1 = np.sqrt(SSB_ab**2 + dd_ab**2 + dw_ab**2 + di_ab**2 + n_ab**2 + n_orb_ab**2 + n_t_ab1**2)
T_ssh_ab2 = np.sqrt(SSB_ab**2 + dd_ab**2 + dw_ab**2 + di_ab**2 + n_ab**2 + n_orb_ab**2 + n_t_ab2**2)
# total relative SSH error
T_ssh_r1 = np.sqrt(SSB_r1**2 + dd_r1**2 + dw_rswot**2 + di_rswot**2 + n_r1**2 + n_orb_r1**2 + n_t_r1**2)
T_ssh_r2 = np.sqrt(SSB_r1**2 + dd_r1**2 + dw_rswot**2 + di_rswot**2 + n_r1**2 + n_orb_r1**2 + n_t_r2**2)
T_ssh_r3 = np.sqrt(SSB_r**2 + dd_r**2 + dw_r**2 + di_r**2 + n_r**2 + n_orb_r**2 + n_t_r**2)

data_screen1 = SSB_data_comb + dd_data_comb + dw_data_comb + di_data_comb + n_data_comb + n_timing_error_comb
data_screen3 = SSB_data_insar + dd_data_insar + dw_data_insar + di_data_insar + n_data_insar + n_timing_error_insar


# PLOT THE RESULT
if np.size(T_ssh_r1)>1:
    plt.figure()
    plt.plot(scale_az/1000, T_ssh_r1)
    plt.plot(scale_az/1000, T_ssh_r3)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Grid spacing in along-track (km)', fontsize=16)
    plt.ylabel('Total accuracy (cm)', fontsize=16)
    plt.title('Total error analysis')
    plt.legend(['Relative error in Concept #1 & #2','Relative error in Concept #3'])
       #    plt.ylim(0, 2)
    plt.grid() 

if ssbplot:
    plt.figure()
    plt.pcolor(gridx/1000, gridy/1000, SSB_data_comb)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Ground-range (km)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('SSB errors for nadir-looking altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)
    plt.figure()
    plt.pcolor(gridx/1000, gridy/1000, SSB_data_insar)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Ground-range (km)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('SSB errors for swath altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)
if ddplot:
    plt.figure()
    plt.pcolor(gridx/1000, gridy/1000, dd_data_comb)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Ground-range (km)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('Dry delay errors for nadir-looking altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)
    plt.figure()
    plt.pcolor(gridx/1000, gridy/1000, dd_data_insar)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Ground-range (km)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('Dry delay errors for swath altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)
if dwplot:
    plt.figure()
    plt.pcolor(gridx/1000, gridy/1000, dw_data_comb)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Ground-range (km)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('Wet delay errors for nadir-looking altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)
    plt.figure()
    plt.pcolor(gridx/1000, gridy/1000, dw_data_insar)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Ground-range (km)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('Wet delay errors for swath altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)
if diplot:
    plt.figure()
    plt.pcolor(gridx/1000, gridy/1000, di_data_comb)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Ground-range (km)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('Ionospheric delay errors for nadir-looking altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)
    plt.figure()
    plt.pcolor(gridx/1000, gridy/1000, di_data_insar)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Ground-range (km)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('Ionospheric delay errors for swath altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)
if alti_noise_plot:
    plt.figure()
    plt.pcolor(gridx/1000, gridy/1000, n_data_comb)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Ground-range (km)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('Altimeter noise for nadir-looking altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)
    plt.figure()
    plt.pcolor(gridx/1000, gridy/1000, n_data_insar)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Ground-range (km)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('Altimeter noise for swath altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)
if orbitplot:
    plt.figure()
    plt.pcolor(np.degrees(grid_swath), gridx/1000, n_orb_data_insar)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Incidence angle (deg)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('Orbit error related error for swath altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)
    plt.figure()
    plt.pcolor(gridx/1000, gridy/1000, n_orb_data_comb)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Ground-range (km)', fontsize=16)
    plt.ylabel('Along-track (km)', fontsize=16)
    plt.title('Orbit error related error for nadir-looking altimeter', fontsize=16)
    plt.grid()
    tt = plt.colorbar()
    tt.set_label('cm', fontsize=16)
    tt.ax.tick_params(labelsize=16)

 