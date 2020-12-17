# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 16:41:51 2019

# THESE FUNCTIONS HELP CALCULATE THE EFFECTIVE LOOK NUMBERS FOR ALTIMETERS

@author: lyh
"""
import numpy as np


def cell_calculate(al_res, raney_l, sar_res):
    if np.size(al_res)>1 and np.size(raney_l)==1:
        cel = np.zeros_like(al_res)
        for i in range(np.size(al_res)):
            if al_res[i]>(2 * raney_l):
                cel[i] = 2 * raney_l / sar_res
            else:
                cel[i] = al_res[i] / sar_res
    elif np.size(raney_l)>1 and np.size(al_res)==1:
        cel = np.zeros_like(raney_l)
        for i in range(np.size(raney_l)):
            if al_res>(2 * raney_l[i]):
                if np.size(sar_res)>1:
                    cel[i] = 2 * raney_l[i] / sar_res[i]
                else:
                    cel[i] = 2 * raney_l[i] / sar_res
            else:
                if np.size(sar_res)>1:
                    cel[i] = al_res / sar_res[i]
                else:
                    cel[i] = al_res / sar_res
    elif np.size(raney_l)>1 and np.size(al_res)>1:
        cel = np.zeros_like(raney_l)
        for i in range(np.size(raney_l)):
            if al_res[i]>(2 * raney_l[i]):
                cel[i] = 2 * raney_l[i] / sar_res[i]
            else:
                cel[i] = al_res[i] / sar_res[i]
    else:
        if al_res>(2 * raney_l):
            cel = 2 * raney_l / sar_res
        else:
            cel = al_res / sar_res
    
    return cel

def area_cal(sar_res, plu_res, pul_rad):
    if np.size(sar_res)>1 and np.size(plu_res)==1:
        A = np.zeros_like(sar_res)
        for i in range(np.size(sar_res)):
            if sar_res[i]>plu_res or sar_res[i]==plu_res:
                A[i] = np.pi * pul_rad**2
            else:
                A[i] = 2 * pul_rad * sar_res[i]
    elif np.size(plu_res)>1 and np.size(sar_res)==1:
        A = np.zeros_like(plu_res)
        for i in range(np.size(plu_res)):
            if sar_res>plu_res[i] or sar_res==plu_res[i]:
                A[i] = np.pi * pul_rad[i]**2
            else:
                A[i] = 2 * pul_rad[i] * sar_res
    elif np.size(plu_res)>1 and np.size(sar_res)>1:
        A = np.zeros_like(sar_res)
        for i in range(np.size(plu_res)):
            if sar_res[i]>plu_res[i] or sar_res[i]==plu_res[i]:
                A[i] = np.pi * pul_rad[i]**2
            else:
                A[i] = 2 * pul_rad[i] * sar_res[i]
    else:
        if sar_res>plu_res or sar_res==plu_res:
            A = np.pi * pul_rad**2
        else:
            A = 2 * pul_rad * sar_res
    return A