#Building Primitive Cells with Hoppings in certain Energy intervals
#Working with the mcell class
#Based on Custom Mos2 Class

import tbplas as tb
import numpy as np
import math
from mos2class import mcell,  metric
import matplotlib.pyplot as plt

amax = 2
bmax = 2
cutoff_distance = 3
#------------------------------------------------------
#1.: Plotting e_min vs remaining hoppings
cell_1 = mcell("cell_1")
allhops_array = cell_1.mhoppingtable()[0]

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
plt.savefig("pca_exp0_nhoppings_newmos2class.png")
#plt.show()
plt.clf()

#------------------------------------------------------
#2.: Adding orbitals with E>Emin

#Experimental cell 2:
cell_2 = mcell("cell_2")
allhops_array, labels = cell_2.mhoppingtable()

#Comparison Cell 3:
cell_3 = mcell("cell_3")

#Loop:
minhopa = 0.0
minhopb = 1.6
minhops = np.linspace(minhopa,minhopb,100)
metrics = []
for minhop in minhops:
    filtered_array = allhops_array[abs(allhops_array[:,6])>minhop]
    cell_2.changehops_toarr(filtered_array)
    metrics.append(metric(cell_2, cell_3))
    print(metric(cell_2,cell_3))

plt.gca().invert_xaxis()
plt.plot(nvec, metrics)
plt.xlabel("Number N of Hoppings / Total Number of Hoppings")
plt.ylabel("Error Metric m")
plt.title("Error Metric m for different N(E_min). \n amax = bmax = %s, cutoff-distance d = %s A" % (amax, cutoff_distance))
plt.savefig("pca_exp0_Metric_vsN_newmos2class.png")
#plt.show()
plt.clf()
