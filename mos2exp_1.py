
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
                plt.plot(k_len[i], bands[i][:, j], color=color, linewidth=width, label = str(expvector[i].name))
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
    plt.title("Band Structure of MoS2 with varying hopping terms. \n cutoff_distance = %s A, minhop = %s eV" % (cutoff_distance, minhop))
    plt.savefig("mos2exp1.png")
    plt.show()
    plt.close()


############
##MAIN:
############

cutoff_distance = 3.5
minhop = 1e-6
expvector = []
for i in range(1,4):
    amax= i
    bmax = i
    tbexp = mos2exp(i)
    tbexp.addhoppings(cutoff_distance,amax,bmax,minhop)
    expvector.append(tbexp)


safebandstructure(expvector)

