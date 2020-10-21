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

E0 = -0.4 # define oxigen onsite energy
# define lattice vectors
lat=[[1.0,0.0],[0.0,1.0]]
r=1.0
t1 = -3.8
t2 = -0.4
# define coordinates of orbitals
# orb=np.zeros((3,2),dtype=float)
orb = [[0, 0] , [1, 0], [0, 1]]
ep = np.zeros(3)
ep[1] = ep[2] = E0

# define model
my_model=tb_model(2,2,lat,orb)
# leave on-site energies to default zero values
# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_onsite([ep[0],ep[1],ep[2]])
my_model.set_hop(t1, 0, 1, [0,0])
my_model.set_hop(t1, 0, 2, [0,0])
my_model.set_hop(t1, 0, 2, [0,-1])
my_model.set_hop(t1, 0, 1, [-1,0])
my_model.set_hop(t2, 1, 2, [0,0])
my_model.set_hop(t2, 1, 2, [1,0])
my_model.set_hop(t2, 1, 2, [0,-1])
my_model.set_hop(t2, 1, 2, [1,-1])
# print tight-binding model
my_model.display()
# define a path in k-space
path = [[0.0, 0.0], [0.5, 0.0], [0.5, 0.5], [0.0, 0.0]]
(k_vec,k_dist,k_node)=my_model.k_path(path,401)
k_label=[r"$G$",r"$X$", r"$M$", r"$G$"]
evals=my_model.solve_all(k_vec)

fig, ax = plt.subplots()
for n in range(3):
    ax.plot(k_dist,evals[n])
ax.set_title("Cu02 band structure")
ax.set_xlabel("K-Path")
ax.set_ylabel("Energy eigenvalue")
ax.set_xticks(k_node)
ax.set_xticklabels(k_label)
# ax.set_xticks(fs_vector)
# ax.set_xlim(fs_vector[0],fs_vector[-1])
fig.tight_layout()
fig.savefig("Cu02_bs.pdf")