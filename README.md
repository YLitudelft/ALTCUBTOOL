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
