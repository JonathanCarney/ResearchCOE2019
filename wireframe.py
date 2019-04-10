# Header
# Author: Jonathan Carney
# Date: 4/9/19
# Title: wireframe.py
# Purpose: produce a wireframe plot of <n> given proposed model for MORF gas absorbtion

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
# import Axes3D as ax

#constants
k = 1.3806*10**-23 #boltzmann's costant
delta = 1 #ring energy
E_heart = 1 #energy as favored site
E_sad = 1 #energy at less favored site

#functions
def beta(T):
    return(1/(k*T))

def microstate_energy(n,r):
    return(r*delta+n*(r*E_heart+(1-r)*E_sad))

def microstate_probability_proportionality(n,r,T,mu):
    return(np.exp(beta(T)*(mu*n-microstate_energy(n,r))))

def partition(T,mu):
    z=0
    for n in range(0,2):
        for r in range(0,2):
            z += microstate_probability_proportionality(n,r,T,mu)
    return(z)

def microstate_probability(n,r,T,mu):
    return(microstate_probability_proportionality(n,r,T,mu)/partition(T,mu))

def expectation_n(T,mu):
    n_exp = 0
    for n in range(0,2):
        for r in range(0,2):
            n_exp += n*microstate_probability(n,r,T,mu)
    return(n_exp)

def expectation_r(T,mu):
    r_exp = 0
    for n in range(0,2):
        for r in range(0,2):
            r_exp += r*microstate_probability(n,r,T,mu)
    return(r_exp)

#producing wirefram plots of <n> and <r>

T_start = 250
T_end = 300
mu_start = 0
mu_end = 100
number_indexes = 100

temp = np.linspace(T_start,T_end,number_indexes)
chem_pot = np.linspace(mu_start,mu_end,number_indexes)

X,Y = np.meshgrid(temp,chem_pot)

#figure 1
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.title("expectation value of n")
ax.plot_wireframe(X,Y,expectation_n(X,Y),color='#b9484e')
plt.xlabel('Temperature(Kelvin)')
plt.ylabel('Chemical Potential(Jule/mol)')
plt.show()

#figure 2
fig2 = plt.figure()
ax = fig2.add_subplot(111, projection='3d')
plt.title("expectation value of r")
ax.plot_wireframe(X,Y,expectation_r(X,Y),color='#647d8e')
plt.xlabel('Temperature(Kelvin)')
plt.ylabel('Chemical Potential(Jule/mol)')
plt.show()

print(expectation_n(X,Y))
