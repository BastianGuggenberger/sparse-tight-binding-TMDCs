
import numpy as np
import tbplas as tb
import matplotlib.pyplot as plt
from numpy.linalg import norm
from math import exp, sqrt

from mos2expclass import mos2exp


def getbandgap(tbexp):

    k_points = np.array([
        [0.0, 0.0, 0.0],    # Gamma
        [1./2, 0.0, 0.0],   # M
        [2./3, 1./3, 0.0],  # K
        [0.0, 0.0, 0.0],    # Gamma
    ])
    k_label = ["G", "M", "K", "G"]
    k_path, k_idx = tb.gen_kpath(k_points, [40, 40, 40])
    #print(len(bands))
    k_len, bands = tbexp.cell.calc_bands(k_path)

    energies = []
    num_k = bands.shape[0]
    num_bands = bands.shape[1]
    
    for i in range(num_bands):
        for j in range(num_k):
            energies.append((bands[j, i]))
    energies = np.array(energies)
    energies = np.sort(energies)
    e_old = energies[0]
    e_new = energies[1]
    deltae = 0
    maxdelta=0
    for i in range(1,len(energies)):
        e_new = energies[i]
        deltae = e_new - e_old
        if(deltae > maxdelta):
            maxdelta = deltae
        e_old = e_new
    return maxdelta

############
##MAIN:
############


amax= 2
bmax = 2
cutoff_distance = 3.5
"""minhop = 1.2

#Bandstructures for different E_min
lowminhop = 1
highminhop = 1.5 #There must be at least 1 Hopping for calcbands to work
minhops = 1.2
tbexp = mos2exp(minhop)
tbexp.addhoppings(cutoff_distance,amax,bmax,minhop)

print("Bandgap = " + str(getbandgap(tbexp)) + " eV")

"""
#Plotting bandwidth vs minhop
lowminhop = 0
highminhop = 2.4
minhops = np.linspace(lowminhop, highminhop,50)
bandgaps = []
for minhop in minhops:
    tbexp = mos2exp(cutoff_distance)
    tbexp.addhoppings(cutoff_distance,amax,bmax,minhop)
    #print(tbexp.hoppings)
    bandgaps.append(getbandgap(tbexp))

plt.plot(minhops, bandgaps,label="Calculated Bandgaps")
plt.axhline(1.77, color='r', linewidth=0.5,label="1.77 eV (Literature)")
plt.xlabel("Energy E_min in eV")
plt.ylabel("Bandgap in eV")
plt.legend()
plt.title("Bandgaps of MoS2 with varying minimal hopping energy. \n amax = bmax = %s, cutoff-distance d = %s A" % (amax, cutoff_distance))
plt.savefig("mos2exp4_bandgaps.png")
#plt.show()
plt.clf()