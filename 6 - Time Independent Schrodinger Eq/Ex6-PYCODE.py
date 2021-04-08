#%% Libraries
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
plt.style.use("default") # Plot white figures


#%% Compile fortran code
    # Checks the existence of the executable, if not it compiles
    # the fortran code
if os.path.isfile('./eigen_ho.x'):
    print('Fortran code is already compiled! \n')
else:
    os.system("gfortran -O2 Ex6-FCODE.f90 -o eigen_ho.x -llapack")
    print("Fortran code has been compiled! \n")

os.system('./eigen_ho.x')


#%% Define theoretical functions for plotting
def eigenvalue(n):
    return (2*n+1)

def QHO(x, number, m=1, omega=1, hbar=1):
    if number==0:
        return (m*omega)**(1/4)*np.exp(-m*omega*x**2/(2*hbar))/(hbar**(1/4)*np.pi**(1/4))
    elif number==1:
        return ((m*omega)/(hbar*np.pi))**(1/4)*np.sqrt(2)*np.sqrt(m*omega/hbar)*x*np.exp(-m*x*x*omega/(2*hbar))
    elif number==2:
        return -((m*omega)/(hbar*np.pi))**(1/4)/np.sqrt(2)*((2*m*x*x*omega/hbar)-1)*np.exp(-m*x*x*omega/(2*hbar))
    elif number==3:
        return -((m*omega)/(hbar*np.pi))**(1/4)/np.sqrt(3)*(2*((m*omega/hbar)**(3/2))*(x**3)-3*np.sqrt(m*omega/hbar)*x)*np.exp(-m*x*x*omega/(2*hbar))
    else:
        return 0


#%% Eigenvalues
ev = np.loadtxt('eigenvalues.txt')
n = np.arange(0,len(ev),1)
th_ev = eigenvalue(n)
print("First Eigenvalues: ", ev[:5])
fig, ax = plt.subplots(1,2, figsize=(20,8))
ax[0].plot(np.arange(0,len(ev)),ev, color='C1', linewidth=5, label='Computed')
ax[0].plot(np.arange(0,len(ev)),th_ev, color='C0', linewidth=5, label='Theoretical')
ax[0].set_title('Eigenvalues')
ax[0].set_xlabel(r'n', fontsize=18)
ax[0].legend(fontsize=18)
ax[1].plot(np.arange(0,30),ev[:30], color='C1', linewidth=5, label='Computed')
ax[1].plot(np.arange(0,30),th_ev[:30], color='C0', linewidth=5, label='Theoretical')
ax[1].set_title('Zoomed Eigenvalues')
ax[1].set_xlabel(r'n', fontsize=18)
ax[1].legend(fontsize=18)
fig.savefig('eigenvalues.png')
plt.show()


#%% Load data and plot
ef = np.loadtxt('eigenvectors.txt')
fig, ax = plt.subplots(len(ef[0,:])-1, 1, figsize=(20,25))
for i in range(0, len(ef[0,:])-1):
    ax[i].plot(ef[:,0], QHO(ef[:,0], i), color=('C'+str(i)), linewidth=5, alpha=0.5, label='Theoretical $\psi_{'+str(i)+'}$')
    ax[i].plot(ef[:,0], ef[:,i+1], color=('C'+str(i)), linestyle='-.', linewidth=5, label=('Computed $\psi_{'+str(i)+'}$'))
    ax[i].set_xlabel(r'x', fontsize=20)
    ax[i].set_ylabel(r'$\psi_{'+str(i)+'}$', fontsize=20)
    ax[i].legend(fontsize=18)
fig.savefig('eigenvectors.png')
plt.show()



# %%
