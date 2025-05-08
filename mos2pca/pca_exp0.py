#Building Primitive Cells with Hoppings in certain Energy intervals
#Working with the mcell class

import tbplas as tb
import numpy as np
import math
from mos2pcaclass import mcell, mhop, metric
import matplotlib.pyplot as plt

amax = 2
bmax = 2
cutoff_distance = 3
#------------------------------------------------------
#1.: Plotting e_min vs remaining hoppings
cell_1 = mcell("cell_1")
allhops_array = cell_1.mfindhoppings(3,2,2,0)[0]

minhopa = 0.0
minhopb = 2.5
total = allhops_array.shape[0]
nvec = []
minhops = np.linspace(minhopa,minhopb,100)

for minhop in minhops:
    n = 0
    for i in range(total):
        if(abs(allhops_array[i][6]) > minhop):
            n += 1
    nvec.append(n/total)
plt.plot(minhops,nvec)
plt.xlabel("Energy E_min in eV")
plt.ylabel("Remaining Hoppings N / Total Hoppings")
plt.title("Hopping Energies of Hopping Terms. \n amax = bmax = %s, cutoff-distance d = %s A" % (amax, cutoff_distance))
plt.savefig("pca_exp0_nhoppings.png")
#plt.show()
plt.clf()

#------------------------------------------------------
#2.: Adding orbitals with E>Emin

#Experimental cell 2:
cell_2 = mcell("cell_2")
allhops_array, labels, allhops_vector = cell_2.mfindhoppings(3,2,2,0)

#Comparison Cell 3:
cell_3 = mcell("cell_3")
cell_3.mnewhoppings(allhops_vector)

#Loop:
metrics = []
for minhop in minhops:
    filtered_array = allhops_array[abs(allhops_array[:,6])>minhop]
    cell_2.changehops_toarr(filtered_array)
    metrics.append(metric(cell_2, cell_3))
    print(metric(cell_2,cell_3))

plt.plot(minhops, metrics)
plt.xlabel("Energy E_min in eV")
plt.ylabel("Error Metric m")
plt.title("Error Metric m for different E_min. \n amax = bmax = %s, cutoff-distance d = %s A" % (amax, cutoff_distance))
plt.savefig("pca_exp0_Metric.png")
#plt.show()
plt.clf()
