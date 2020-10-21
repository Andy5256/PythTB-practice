#!/usr/bin/env python

# Tight-binding model for p_z states of benzene molecule with 
# uniform electric field along x, and the site energies be raised
# or lowered in a way that is linear in their spatial coordinates.

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)

from __future__ import print_function # python3 style print
from pythtb import * # import TB model class

import numpy as np
import matplotlib.pyplot as plt

# define lattice vectors
lat=[[1.0,0.0],[0.0,1.0]]
r=1.2
ep = np.zeros(6)
fs_vector = []
eval_vector = []
for fs in range(10) :
    fs = fs*0.2 # electric field strenth
# define coordinates of orbitals
    orb=np.zeros((6,2),dtype=float)
    for i in range(6) :
        angle=i*np.pi/3.0
        orb[i,:] = [r*np.cos(angle) , r*np.sin(angle)]
    for j in range(6) :
        angle=j*np.pi/3.0
        ep[j] = -0.4 + fs*np.cos(angle)
    t=-0.25
# define model
    my_model=tb_model(0,2,lat,orb)
# leave on-site energies to default zero values
# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
    my_model.set_onsite([ep[0],ep[1],ep[2],ep[3],ep[4],ep[5]])
    my_model.set_hop(t,0,1)
    my_model.set_hop(t,1,2)
    my_model.set_hop(t,2,3)
    my_model.set_hop(t,3,4)
    my_model.set_hop(t,4,5)
    my_model.set_hop(t,5,0)
# print tight-binding model
    my_model.display()

    (eval,evec)=my_model.solve_all(eig_vectors=True)
    fs_vector.append(fs)
    for k in range(6) :
        eval_vector.append(eval[k])
# print (" %i" % len(eval_vector))
a = np.mat(eval_vector)
b = a.reshape(10,6)
# print (type(b))
    # np.set_printoptions(formatter={'float':'{: 6.3f}'.format})
    # print(" n eigval eigvec")
    # for n in range(6) :
    #     print (" %2i %7.3f  " % (n,eval[n]), evec[n,:].real)
    # plot band structure
fig, ax = plt.subplots()
for n in range(6):
    ax.plot(fs_vector,b[:,n])
ax.set_title("Eval with respect to field strength")
ax.set_xlabel("Electric field strenth")
ax.set_ylabel("Energy eigenvalue")
ax.set_xticks(fs_vector)
ax.set_xlim(fs_vector[0],fs_vector[-1])
fig.tight_layout()
fig.savefig("fs_eval.pdf")