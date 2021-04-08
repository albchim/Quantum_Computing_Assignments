#%% Libraries
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
plt.style.use("default") # Plot white figures

#%% Compile fortran code
    # Checks the existence of the executable, if not it compiles
    # the fortran code
if os.path.isfile('./ho_ev.x'):
    print('Fortran code is already compiled! \n')
    decision=(input("Recompile fortran script? [y/n]")).lower()
    if decision=='y':
        os.system("gfortran -O2 Ex7-FCODE.f90 -llapack -lfftw3 -o ho_ev.x")
        print("Fortran code has been compiled! \n")
    if decision!='y' and decision!='n':
        print('Invalid input')
        exit()
    if decision=='n':
        print('Fortran code is already compiled! \n')
else:
    os.system("gfortran -O2 Ex7-FCODE.f90 -llapack -lfftw3 -o ho_ev.x")
    print("Fortran code has been compiled! \n")


#%% Generates input file and runs code

decision=(input("Run fortran script? [y/n]")).lower()
if decision=='y':
    dim = open('input.txt', 'w')
    n_time=int(input("Enter the number of time steps to perform (INTEGER NUMBER): "))
    t_max=int(input("Enter the value of T_max (INTEGER NUMBER): "))
    dim.write(str(n_time)+"\n")
    dim.write(str(t_max)+".\n")
    dim.close()
    os.system('./ho_ev.x < input.txt')

if decision!='y' and decision!='n':
    print('Invalid input')
    exit()
if decision=='n':
    print('Using already existing data')

#%% Load data inputs
dim = open('input.txt', 'r')
inputs = (dim.readlines())
n_time = int(inputs[0])
t_max = int(inputs[1][:len(inputs[1])-2])


#%% Load data
from scipy.integrate import simps

t_ev = np.loadtxt('t_evolution.txt')
print(t_ev.shape)
fig, ax = plt.subplots(figsize=(10,5))
ax.set_xlim([-3,3])
n_curves = 20
i_max = max(len(t_ev[0,:])-1, n_curves)
if i_max==len(t_ev[0,:])-1:
    i_step = int((len(t_ev[0,:])-1)/n_curves)
else:
    i_step = 1
for i in np.arange(0, len(t_ev[0,:])-1, i_step):
    if i%100==0:
        print("Integral of |\psi_{0,"+str(i)+"}|^2:",simps(t_ev[:,i+1], t_ev[:,0]))
        ax.plot(t_ev[:,0], t_ev[:,i+1], linestyle='-', linewidth=2, label=('$\psi_{0,'+str(i)+'}$'))
    elif i%100!=0:
        ax.plot(t_ev[:,0], t_ev[:,i+1], linestyle='-', linewidth=2)
ax.plot(t_ev[:,0], t_ev[:,(len(t_ev[0,:]))-1], color="palevioletred", linestyle='-', linewidth=3, label=('$\psi_{0,'+str((len(t_ev[0,:])-1))+'}$'))
ax.set_xlabel(r'x', fontsize=20)
ax.set_title(r"dt ="+str(t_max/n_time)+r'    $T_{max}=$'+str(t_max), fontsize=20)
ax.set_ylabel(r'$\psi_{0}$', fontsize=20)
ax.legend()#fontsize=18)
fig.savefig('t_evolution.png')
plt.show()
plt.close()


#%% Gif plotting
decision = 'y'
print("\n Make a gif plot?\n WARNING: there may be some difficulties with libraries and visualization....")
decision=(input("Type [y/n]")).lower()


if decision=='y':
    import matplotlib.animation as animation

    fig, ax = plt.subplots(figsize=(10,5))

    x = t_ev[:,0]
    line, = ax.plot(x, t_ev[:,i], "C0", linewidth=2)
    ax.set_ylim([0,0.8])
    ax.set_xlim([-4,5])

    def animate(i):
        line.set_ydata(t_ev[:,i])  # update the data
        return line,

    # Init only required for blitting to give a clean slate.
    def init():
        line.set_ydata(t_ev[:,1])
        return line,

    ani = animation.FuncAnimation(fig, animate, frames=np.arange(1, len(t_ev[0,:])-1, i_step), init_func=init, interval=200, blit=True)
    ani.save("t_evolution.gif")
else:
    exit()

