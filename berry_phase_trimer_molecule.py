#!/usr/bin/env python

# Tight-binding model for 3.4 Berry phase of the ground state of the trimer molecule.

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)

from __future__ import print_function # python3 style print
from pythtb import * # import TB model class

import numpy as np
import matplotlib.pyplot as plt

# define lattice vectors
lat=[[1.0,0.0],[0.0,1.0]]
r=1
evec_vector = []
eval_vector = []
berry_phase = 0
print(berry_phase)
# define coordinates of orbitals
orb=np.zeros((3,2),dtype=float)
for i in range(3) :
    # evec_vector = []
    # eval_vector = []
    angle=i*2*np.pi/3.0+np.pi/2.0
    # print(angle)
    orb[i,:] = [r*np.cos(angle) , r*np.sin(angle)]
# print(orb)
# define model
# for fi in range(1) :
# for fi in 0, 1, 2 :
for fi in 0, 2*np.pi/3.0, 4*np.pi/3.0 :
    my_model=tb_model(0,2,lat,orb)
# leave on-site energies to default zero values
# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
    my_model.set_onsite([0.0, 0.0, 0.0])   
    t0=-1.0
    s=-0.4
    t01 = t0 + s*np.cos(fi)
    t12 = t0 + s*np.cos(fi-2*np.pi/3.0)
    t20 = t0 + s*np.cos(fi-4*np.pi/3.0)
    alpha = -np.pi/4
    # my_model.set_hop(t01, 0, 1)
    # my_model.set_hop(t12, 1, 2)
    # my_model.set_hop(t20, 2, 0)
    my_model.set_hop(t01, 0, 1)
    my_model.set_hop(t12, 1, 2)
    my_model.set_hop(t20*np.exp((0+1j)*alpha), 2, 0)
    # print(t20*np.exp((0+1j)*alpha/3))
# print tight-binding model
    # my_model.display()
    (eval,evec)=my_model.solve_all(eig_vectors=True)
    print(evec[0])
    for k in range(1) :
        evec_vector.append(evec[k])
        eval_vector.append(eval[k])
# print(evec_vector)
    # print(eval_vector)
# print (" %i" % len(eval_vector))
a = np.mat(evec_vector)
b = a.reshape(3,3)
print(berry_phase)
# print(b)
# b = np.transpose(b)
print(b)
b0 = b[0]
b1= b[1]
b2 = b[2]
print(b0)
print(b1)
print(b2)
prod = np.vdot(b0,b1)*np.vdot(b1,b2)*np.vdot(b2,b0)
# print(prod)
# berry_phase = -np.angle(prod)
berry_phase = -np.angle(np.vdot(b0,b1)) -np.angle(np.vdot(b1,b2)) -np.angle(np.vdot(b2,b0))
print(np.angle(np.vdot(b1,b2)))
print(berry_phase)
    # np.set_printoptions(formatter={'float':'{: 6.3f}'.format})
    # print(" n eigval eigvec")
    # for n in range(6) :
    #     print (" %2i %7.3f  " % (n,eval[n]), evec[n,:].real)
    # plot band structure
# fig, ax = plt.subplots()
# for n in range(6):
#     ax.plot(delta_vector,b[:,n])
# ax.set_title("Eval with respect to delta")
# ax.set_xlabel("\delta")
# ax.set_ylabel("Energy eigenvalue")
# ax.set_xticks(delta_vector)
# ax.set_xlim(delta_vector[0],delta_vector[-1])
# fig.tight_layout()
# fig.savefig("delta_eval.pdf")