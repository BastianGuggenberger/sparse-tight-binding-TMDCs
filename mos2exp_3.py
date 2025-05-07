
import numpy as np
import tbplas as tb
import matplotlib.pyplot as plt
from numpy.linalg import norm
from math import exp, sqrt

from mos2expclass import mos2exp


def safebandstructure(expvector,lowminhop,highminhop):
    n= len(expvector)
    cellvector = [exp.cell for exp in expvector]

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
    #print(len(bands))
    for i in range(n):
        k_lenn, bandss = cellvector[i].calc_bands(k_path)
        k_len.append(k_lenn)
        bands.append(bandss)

    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628"]

    #Plotting the Bandstructure
    for i in range(n):
        color = colors[i]
        num_bands = bands[i].shape[1]
        width = 0.7
        for j in range(num_bands):
            if(j==0):
                label = str(expvector[i].name) + "eV, " + str(expvector[i].hoppings) +" Hoppings"
                plt.plot(k_len[i], bands[i][:, j], color=color, linewidth=width, label = label)
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
    plt.title("Band Structure of MoS2 for different E_min\n amax = bmax = %s, Cutoff-Distance d = %s A" % (amax, cutoff_distance))
    pngname = "mos2exp_3_l_low" + str(lowminhop) + "_high" + str(highminhop) + ".png"
    plt.savefig(pngname)
    #plt.show()
    plt.close()




############
##MAIN:
############


amax= 2
bmax = 2
cutoff_distance = 3.5

#Plotting n_hoppings vs minhop
lowminhop = 0
highminhop = 2.7
minhops = np.linspace(lowminhop, highminhop,50)
n_hoppings = []
for minhop in minhops:
    tbexp = mos2exp(cutoff_distance)
    tbexp.addhoppings(cutoff_distance,amax,bmax,minhop)
    #print(tbexp.hoppings)
    n_hoppings.append(tbexp.hoppings)

plt.plot(minhops, n_hoppings)
plt.xlabel("Energy E_min in eV")
plt.ylabel("Number of Hoppings N with E > E_min")
plt.title("Number of Hoppings of MoS2 with varying minimal hopping energy. \n amax = bmax = %s, cutoff-distance d = %s A" % (amax, cutoff_distance))
plt.savefig("mos2exp3_nhoppings.png")
#plt.show()
plt.clf()



#Bandstructures for different E_min
lowminhop = 0
highminhop = 0.05 #There must be at least 1 Hopping for calcbands to work
minhops = np.linspace(lowminhop, highminhop,5)
expvector = []

for minhop in minhops:
    tbexp = mos2exp(minhop)
    tbexp.addhoppings(cutoff_distance,amax,bmax,minhop)
    print(tbexp.hoppings)
    expvector.append(tbexp)

safebandstructure(expvector,lowminhop,highminhop)