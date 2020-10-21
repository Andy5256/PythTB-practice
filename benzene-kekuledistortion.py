#!/usr/bin/env python

# Tight-binding model for p_z states of benzene molecule with 
# kekule distortion.

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)

from __future__ import print_function # python3 style print
from pythtb import * # import TB model class

import numpy as np
import matplotlib.pyplot as plt

# define lattice vectors
lat=[[1.0,0.0],[0.0,1.0]]
r=1.2
delta_vector = []
eval_vector = []
for delta in range(10) :
    delta = delta*0.05 # electric field strenth
# define coordinates of orbitals
    orb=np.zeros((6,2),dtype=float)
    for i in range(6) :
        angle=i*np.pi/3.0
        orb[i,:] = [r*np.cos(angle) , r*np.sin(angle)]
    ep = -0.4
    t=-0.25
# define model
    my_model=tb_model(0,2,lat,orb)
# leave on-site energies to default zero values
# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
    my_model.set_onsite([ep,ep,ep,ep,ep,ep])
    my_model.set_hop(t+delta,0,1)
    my_model.set_hop(t-delta,1,2)
    my_model.set_hop(t+delta,2,3)
    my_model.set_hop(t-delta,3,4)
    my_model.set_hop(t+delta,4,5)
    my_model.set_hop(t-delta,5,0)
# print tight-binding model
    my_model.display()

    (eval,evec)=my_model.solve_all(eig_vectors=True)
    delta_vector.append(delta)
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
    ax.plot(delta_vector,b[:,n])
ax.set_title("Eval with respect to delta")
ax.set_xlabel("\delta")
ax.set_ylabel("Energy eigenvalue")
ax.set_xticks(delta_vector)
ax.set_xlim(delta_vector[0],delta_vector[-1])
fig.tight_layout()
fig.savefig("delta_eval.pdf")