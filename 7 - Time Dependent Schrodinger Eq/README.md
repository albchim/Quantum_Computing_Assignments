Needed python libraries: numpy, scipy, matplotlib, sys, os
Needed fortran libraries: LAPACK, FFTW3

Simply run the python script with the command:

 >python3 Ex7-PYCODE.py

The user will be prompted to the choice of recompiling or running the fortran executable 
and generate the data files needed for the graphical analysis.
It will also propt the choice for an optional animated gif plot to be produce as a better
data visualization method.


In order to reproduce the results present in the report, just move the 'input.txt' file and the correspondingly named .txt file
in the main code folder and rename the latter as 't_evolution.txt'. Then run the python script and select
'n' when prompted to the choice of running the fortran script so that it can use pre existing data. 