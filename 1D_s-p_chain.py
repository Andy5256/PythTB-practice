#!/usr/bin/env python

# Tight-binding model for 3.4 Berry phase of the ground state of the trimer molecule.

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)

from __future__ import print_function # python3 style print
from pythtb import * # import TB model class

import numpy as np
import matplotlib.pyplot as plt

# define lattice vectors
lat=[[1.0]]
# define coordinates of orbitals
orb=[[0.0],[0.0]]
# print(orb)
# define model
my_model=tb_model(1,1,lat,orb)
# leave on-site energies to default zero values
# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_onsite([0.0, 8.0])   
for t in 0.5, 0.8, 1.0, 1.2 :
    Vss=-1.40*t
    Vpp=3.24*t
    Vsp=1.84*t
    my_model.set_hop(Vss, 0, 0, [1])
    my_model.set_hop(Vsp, 0, 1, [1])
    my_model.set_hop(Vpp, 1, 1, [1])
    my_model.set_hop(-Vsp, 1, 0, [1])
# define a path in k-space
    path = [[0.0], [0.5]]
    (k_vec,k_dist,k_node)=my_model.k_path(path,401)
    k_label=[r"$\gamma$",r"$X$"]
# print tight-binding model
    my_model.display()
    evals=my_model.solve_all(k_vec)
# plot band structure
    fig, ax = plt.subplots()
    ax.plot(k_dist,evals[0])
    ax.set_title("1D chain band structure")
    ax.set_xlabel("Path in k-space")
    ax.set_ylabel("Band energy")
    ax.set_xticks(k_node)
    ax.set_xticklabels(k_label)
    ax.set_xlim(k_node[0],k_node[-1])
    for n in range(len(k_node)):
        ax.axvline(x=k_node[n], linewidth=0.5, color='k')
    fig.tight_layout()
    fig.savefig("1D_band.pdf")