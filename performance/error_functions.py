# -*- coding: utf-8 -*-
"""
Created on Mon 30 10:42:36 2020
It is the file generating error sources for Alticubes
@author: lyh
"""
# Setting up the packages
import numpy as np
from matplotlib import pyplot as plt
from constants import constants as const

# Exponential distribution for altimeter noise
def altimeter_rand(num_x, num_y, sigma):
    speckle = (1./np.sqrt(2) * (np.random.normal(0., 1., size=[num_x, num_y]) +
               1j * np.random.normal(0., 1., size=[num_x, num_y]))).astype(np.complex64)
    amp = np.abs(speckle)
    intensity = amp**2
    intensity = sigma / np.std(intensity) * intensity 
    return intensity

# Interferogram speckle for swath altimeter
def altimeter_rand_insar(num_x, num_y):
    speckle1 = (1./np.sqrt(2) * (np.random.normal(0., 1., size=[num_x, num_y]) +
               1j * np.random.normal(0., 1., size=[num_x, num_y]))).astype(np.complex64)
    speckle2 = (1./np.sqrt(2) * (np.random.normal(0., 1., size=[num_x, num_y]) +
               1j * np.random.normal(0., 1., size=[num_x, num_y]))).astype(np.complex64)
    insar_phase = np.angle(speckle1*np.conj(speckle2))
    return insar_phase

# Error sources with power-law behavior and relative error
def re_acuracy(dis, sig, index, Nx=1000, dx=1):
    kx =  np.fft.fftfreq(Nx, dx)    
#    f1 = (np.abs(kx))>1e-7
    dkx = kx[1] - kx[0]
    spc = ((kx**2)**(index/2)) * dkx
    spc[0] = spc[1]
#    spc = spc * f1
    factor =  sig**2 / np.sum(spc)
    spc = spc * factor
    R = np.abs(np.fft.ifft(spc))
    t = R[int(dis / dx)] / R.max()
    vara = (2 - 2 * t) * sig**2
    return vara

# Relative error of residual wet delay when 1-beam radiometer is applied
def re_wet_one_radiometer(dis, sig, p1, p2, Nx=1000, dx=1):
    kx =  np.fft.fftfreq(Nx, dx)   
    dkx = kx[1] - kx[0]
    spc1 = (9.5e-5 * (kx**2)**(p1/2)) * dkx * (np.abs(kx)>=1e-3) * (np.abs(kx)<0.0023)
    spc2 = (0.036 * (kx**2)**(p2/2)) * dkx * (np.abs(kx)>=0.0023)* (np.abs(kx)<0.0683)
    spc3 = 0.32*(np.abs(kx)>=0.0683) *dkx
    spc = spc1 + spc2 + spc3
    spc = spc - 0.32 * dkx
    spc[0] = spc[1]
#    spc = spc * f1
    factor =  sig**2 / np.sum(spc)
    spc = spc * factor
    R = np.abs(np.fft.ifft(spc))
    t = R[int(dis / dx)] / R.max()
    vara = (2 - 2 * t) * sig**2
    return vara

# Generating the error data with a power-law spectrum
def re_acuracy_data(sig, p=-8/3, Nx=4001, Ny=4001, thr=1e-7, dx=1100):
    kx = np.fft.fftfreq(Nx, dx) + 1e-50
    ky = np.fft.fftfreq(Ny, dx) + 1e-50
    f1 = (np.abs(kx))>thr
    f2 = (np.abs(ky))>thr
    f11, f22 = np.meshgrid(f1,f2)
    dkx = kx[1] - kx[0]
    dky = ky[1] - ky[0]
    kx = kx.reshape((Nx, 1))
    ky = ky.reshape((1, Ny))
    spc = ((kx**2 + ky**2)**(p/2)) * dkx * dky
    f11 = np.fft.fftshift(f11)
    f22 = np.fft.fftshift(f22)
    spc = np.fft.fftshift(spc)
    spc = spc * np.transpose(f11 *f22)
    spc = np.fft.fftshift(spc)
    factor =  sig**2 / np.sum(spc)
    spc = spc * factor
    wn = np.random.randn(Nx, Ny) + 1j * np.random.randn(Nx, Ny)
    data = wn * Nx * Ny * np.sqrt(spc)
    data = np.fft.ifft2(data).real
    return data

# Error data of the orbit error
def re_acuracy_orbit_data(std, gridx, gridy, scalex=40000e3, scaley=10):
    R = std**2 * np.exp(-np.abs(gridx.reshape(gridx.size,1) / scalex) - np.abs(gridy.reshape(1, gridy.size) / scaley))
    P = np.abs(np.fft.fft2(R)) / gridx.size / gridy.size
#    A = std**2 / (1/(res*num)) * (1/(resa*num)) / (np.sum(np.sum(P)))
#    P = P * A
    wn = np.random.randn(gridx.size, gridy.size) + 1j * np.random.randn(gridx.size, gridy.size)
    screen_k = wn * gridx.size * gridy.size * np.sqrt(P)
    data = np.fft.ifft2(screen_k).real
    kk = std/np.std(data)
    data = data * kk        
    return data

# Error data of the orbit slope
def re_accuracy_orbit_slope(grid, grid_swath, kf):
    s = np.ones_like(grid)
    factor = grid_swath * kf
    data = s.reshape((s.size,1)) * factor
    return data

# Random timing error
def timing_error_data(num_x, num_y, sigma):
    data = np.random.normal(0., sigma, size=[num_x, num_y])
    return data    