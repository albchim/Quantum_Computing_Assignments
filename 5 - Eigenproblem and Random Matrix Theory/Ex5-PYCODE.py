# -*- coding: utf-8 -*-

#%% Importing libraries
import numpy as np
import subprocess
import sys
import os

#%% Compile fortran code
    # Checks the existence of the executable, if not it compiles
    # the fortran code
if os.path.isfile('./eigen_spacing.x'):
    print('Fortran code is already compiled! \n')
else:
    os.system("gfortran -O2 Ex5-FCODE.f90 -o eigen_spacing.x -llapack")
    print("Fortran code has been compiled! \n")


#%% Generates input

decision=(input("Run fortran script? [y/n]")).lower()
if decision=='y':
    dim = open('input.txt', 'w')
    dimension=(input("Insert desired matrix dimension to run the analysis: "))
    nruns=int(input("Insert the desired number of sampling runs: "))
    dim.write(dimension)
    dim.close()
    # Create output files from scratch/overwrite the existing ones
    out_h = open("out_h.txt","w")
    out_d = open("out_d.txt","w")
    for i in range(nruns):
        print("----------> Run: ", i)
        # Loops over the number of runs desired and run the fortran script
        os.system('./eigen_spacing.x < input.txt')
    out_h.close()
    out_d.close()

if decision!='y' and decision!='n':
    print('Invalid input')
    exit()
if decision=='n':
    print('Using already existing data')


#%% Plotting
os.system('gnuplot Ex5-GNUCODE.gnu')