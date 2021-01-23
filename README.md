# ALTCUBTOOL
It is a tool to conduct the system requirement and performance analyses of radar altimeters

## Functions & Purpose
It provides:

Code to calculate the required power, antenna size, and data rate for radar altimeters.
Code to conduct an analysis of SSH product performances under various error sources including sea state bias, atmospheric delay, orbit error, altimeter noise, and timing error.
Code to simulate different altimeter constellation concepts: Comb, Specular, and Interferometric swath constellations.
The provided code is flexible and will be frequently updated. They can also be adapted to the analyses for other spaceborne radar altimeter concepts. 

## A bit of history

Some requirement analysis part is from the basic algorithms in an excel file programmed by Peter Hoogeboom (GRS TU Delft). 

In 2019, the code becomes a basic tool to conduct the requirement and performance analysis for a European Space Agency (ESA) project called ALTICUBES by TU Delft.

## General Roadmap

It begins with the basic requirement analysis for power budget, data rate, and other system parameters. Then, it goes to the SSH performance analysis, which is the most important product for radar altimeters. It also includes the performance analysis for the generated swaths.

## Things to be done

Some other codes for calculating ground tracks of Cubesat altimeter constellations cooperating with STK software are in Matlab code, which can be extended in a python version in the future.

## Usage of the files

# Performance package
1) SSH_accuracy_performance.py
It is used to calculate the absolute SSH and relative SSH accuracy. Three constellations are calculated as explained in the heading of the file.
The error spectra considered here are mainly from the SWOT system analysis. The reference is: [RD].	Daniel Esteban Fernandez, et al. SWOT Project Mission Performance and Error Budget, JPL D-79084, 2017. 

Input: The input parameters are in the cfg file (ssh part). If you want to change the values of the absolute SSH accuracy after different correction methods, you can do it in the cfg file.
Output: 
Calculated absolute accuracy (T_ssh_ab1 and T_ssh_ab2) 
Relative accuracy (T_ssh_r1, T_ssh_r2 and T_ssh_r3). 
1: comb constellation; 2: specular constellation; 3: interferometric constellation.

2) plot_scattering.py
It is used to generate the curve between NRCS and the wind speed. The NRCS data is simulated by OceanSAR simulator. You can directly run it to plot the data. All the scattering data is provided in the folder of scattering.

Input: NRCS data in the scattering folder
Output: the curve between NRCS and the wind speed.

3) Swath_altimeter.py
It conducts the swath performance analysis in the interferometric constellation.

Input: when the off_nadir is a vector
Output: relationship between the off nadir angle and the swath length.

Input: when alt and off nadir are vectors
Output: relationship between the orbit height and the swath length under different off nadir angles

Input: when alt is a vector and off nadir not
Output: relationship between the cross track antenna size of the Microsat and the swath length under different orbit height.

4) height_ambiguity.py
It analyzes the height ambguity under different baseline lengths and incidence angles.

Input: when the incidence angle is a vector
Output: relationship between the baseline length and the height ambguity under different incidence angles.

Input: when the orbit height is a vector
Output: relationship between the baseline length and the height ambguity under different orbit heights.

5) error_function.py
It is a error function paackage to generate different error components in SSH measurements based on the SWOT error spectra. The reference is provided in [RD].
You can change the error spectra if you prefer other error spectra in SSH measurements.

6) Critical_baseline_calculation.py
This file has multiple functions.
First, it is used to calculate the critical baseline in the interferometric constellation.
Second, it can be used to calculate the cross track resolution and baseline roll angle error in the interferometric constellation.
Third, it calculate the gerometric performance (distance between  cubesats).

Input: when the microsat antenna length and the bandwidth are vectors
Output: relationship between the bandwidth and the max. critical baseline under different microsat antenna length.

Input: when the incidence angle and the bandwidth are vectors
Output: relationship between the bandwidth and the max. critical baseline under different incidence angles.

Input: when the incidence angle is a vector
Output: relationship between the incidence angle and the altimeter noise.

Input: when the incidence angle is a vector
Output: the distance between the microsat and cubesat as a function of the incidence angle at near, far, and average range.

Input: when the incidence angle and the bandwidth are vectors
Output: relationship between the incidence angle and the baseline roll angle error under different bandwidths.





