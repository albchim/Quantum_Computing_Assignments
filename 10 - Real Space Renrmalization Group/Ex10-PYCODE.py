#%% Libraries
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
plt.style.use("default") # Plot white figures


#%% Compile fortran code
    # Checks the existence of the executable, if not it compiles
    # the fortran code
if os.path.isfile('./ising-renorm.x'):
    print('Fortran code is already compiled! \n')
    decision=(input("Recompile fortran script? [y/n]")).lower()
    if decision=='y':
        os.system("gfortran -O2 Ex10-FCODE.f90 -o ising-renorm.x -llapack")
        print("Fortran code has been compiled! \n")
    if decision!='y' and decision!='n':
        print('Invalid input')
        exit()
    if decision=='n':
        print('Fortran code is already compiled! \n')
else:
    os.system("gfortran -O2 Ex10-FCODE.f90 -o ising-renorm.x -llapack")
    print("Fortran code has been compiled! \n")

#%% Generates input file and runs code

decision=(input("Run fortran script? [y/n]")).lower()
if decision=='y':
    dim = open('input.txt', 'w')
    N=int(input("Enter the number of qbits for each subsystem (INTEGER NUMBER): "))
    renormiter=int(input("Enter the number of desired subsystems to add (INTEGER NUMBER): "))
    dim.write(str(N)+"\n")
    dim.write(str(renormiter)+"\n")
    dim.close()
    os.system('./ising-renorm.x < input.txt')

if decision!='y' and decision!='n':
    print('Invalid input')
    exit()
if decision=='n':
    print('Using already existing data')


#%% Load data inputs
dim = open('input.txt', 'r')
inputs = (dim.readlines())
N = int(inputs[0])
d = int(inputs[1])


#%% Load data and plot
ef = np.loadtxt('l_eigenvalues.txt')
fig, ax = plt.subplots(figsize=(10,5))
for i in range(0, len(ef[0,:])-1):
    ax.plot(ef[:,0], ef[:,i+1]/N**(d+1), color=('C'+str(i)), linestyle='-', marker='.', linewidth=2, label=('$E_{'+str(i)+'}$'))
    #ax.scatter(ef[:,0], ef[:,i+1], '.', s=3)
    ax.set_xlabel(r'$\lambda$', fontsize=18)
    ax.set_ylabel(r'E($\lambda$)', fontsize=18)
    ax.axvline(0,color='C3',linewidth=2)
    ax.legend(fontsize=12)
    ax.set_title('N='+str(N)+" d="+str(d), fontsize=12)
    ax.set_xlim(min(ef[:,0])-0.5, max(ef[:,0])+0.5)
fig.savefig(('N='+str(N)+" d="+str(d)+'_gstate_energy.png'))
plt.show()