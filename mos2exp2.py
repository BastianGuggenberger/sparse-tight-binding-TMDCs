
import numpy as np
import tbplas as tb
import matplotlib.pyplot as plt
from numpy.linalg import norm
from math import exp, sqrt

from mos2expclass import mos2exp


def safebandstructure(expvector):
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
        width = 1.5
        for j in range(num_bands):
            if(j==0):
                label = str(expvector[i].name) + "A, " + str(expvector[i].hoppings) +" Hoppings"
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
    plt.title("Band Structure of MoS2 with varying cuttoff_distance (in A). \n amax = bmax = %s, minhop = %s eV" % (amax, minhop))
    plt.savefig("mos2exp2.png")
    #plt.show()
    plt.close()




############
##MAIN:
############

minhop = 1e-5
amax= 2
bmax = 2


#Plotting cutoff_distance vs n_hoppings
lowcutoff = 0.1
highcutoff = 1.5
cutoffs = np.linspace(lowcutoff, highcutoff,50)
n_hoppings = []
for cutoff_distance in cutoffs:
    tbexp = mos2exp(cutoff_distance)
    tbexp.addhoppings(cutoff_distance,amax,bmax,minhop)
    n_hoppings.append(tbexp.hoppings)

plt.plot(cutoffs, n_hoppings)
plt.xlabel("Cutoff Distance d (A)")
plt.ylabel("Number of Hoppings")
plt.title("Number of Hoppings of MoS2 with varying cuttoff_distance (in A). \n amax = bmax = %s, minhop = %s eV" % (amax, minhop))
plt.savefig("mos2exp2_nhoppings.png")
plt.clf()
#plt.show()


expvector = []
lowcutoff = 0.5
highcutoff = 1.3
cutoffs = [0.4,0.6,0.8,1,1.2] #With these Distances we cover all the variability in n_hopping
for cutoff_distance in cutoffs:
    tbexp = mos2exp(cutoff_distance)
    tbexp.addhoppings(cutoff_distance,amax,bmax,minhop)
    expvector.append(tbexp)

safebandstructure(expvector)
