#Calculation of Graphene Band Structure

import tbplas as tb
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from numpy.linalg import norm


#################
#PRIMITIVE CELLS:
#Building 7 Primitive Cells, any of these will receive different hoppings later.

a = 2.46
graphene_vectors = np.array([
    [a,0.000,0.000],
    [a*0.5,a*0.5*math.sqrt(3),0.000],
    [0.000,	0.000,	10.0]
])

atom1loc = 1/3 * graphene_vectors[0] +1/3 * graphene_vectors[1]
atom2loc = 2/3 * graphene_vectors[0] +2/3 * graphene_vectors[1]
graphene_atom_locations = np.array([atom1loc,atom2loc])

graphene_cell = []
combinations = [[0],[1],[2],[0,1],[0,2],[1,2],[0,1,2]]
n = len(combinations)
for i in range(n):
    graphene_cell.append(tb.PrimitiveCell(graphene_vectors))


#########
#ORBITALS
#All Cells have the two pz orbitals of monolayer graphene:
for i in range(n):
    graphene_cell[i].add_orbital_cart(graphene_atom_locations[0], unit=tb.ANG, energy=0.0, label="pz")
    graphene_cell[i].add_orbital_cart(graphene_atom_locations[1], unit=tb.ANG, energy=0.0, label="pz")




#Hoppings:
def addhopping(cell,i):
    print(i)
    if(i==0):
        cell.add_hopping(rn=[0, 0], orb_i=0, orb_j=1, energy=-2.7)
    if(i==1):
        cell.add_hopping(rn=[1, 0], orb_i=1, orb_j=0, energy=-2.7)
    if(i==2):
        cell.add_hopping(rn=[0, 1], orb_i=1, orb_j=0, energy=-2.7)

for i in range(n):
    for j in range(len(combinations[i])):
        addhopping(graphene_cell[i],combinations[i][j])


for i in range(n):
    graphene_cell[i].plot()



k_points = np.array([
    [0.0, 0.0, 0.0],    # Gamma
    [1./2, 0.0, 0.0],   # M
    [2./3, 1./3, 0.0],  # K
    [0.0, 0.0, 0.0],    # Gamma
])
k_label = ["G", "M", "K", "G"]
k_path, k_idx = tb.gen_kpath(k_points, [40, 40, 40])
k_len=[]
bands=[]
print(len(bands))
for i in range(n):
    k_lenn, bandss = graphene_cell[i].calc_bands(k_path)
    k_len.append(k_lenn)
    bands.append(bandss)

colors =["yellow","orange","brown","cyan","blue", "green", "red"]

#Plotting the Bandstructure
for i in range(n):
    color = colors[i]
    num_bands = bands[i].shape[1]
    if(i==6):
        width = 2.0
    else:
        width = 1.0
    for j in range(num_bands):
        if(j==0):
            plt.plot(k_len[i], bands[i][:, j], color=color, linewidth=width, label = str(combinations[i]))
        else:
            plt.plot(k_len[i], bands[i][:, j], color=color, linewidth=width)

    for idx in k_idx:
        plt.axvline(k_len[i][idx], color='k', linewidth=1.0)
    plt.xlim((0, np.amax(k_len[i])))
    plt.xticks(k_len[i][k_idx], k_label)
plt.xlabel("k (1/nm)")
plt.ylabel("Energy (eV)")
#plt.tight_layout()
plt.legend()
plt.title("Band Structure of ML-Graphene with varying hopping terms")
plt.show()
plt.close()


