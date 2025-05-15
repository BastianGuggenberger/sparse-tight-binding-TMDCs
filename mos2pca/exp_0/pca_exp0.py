#Building Primitive Cells with Hoppings in certain Energy intervals
#Working with the mcell class
#Based on Custom Mos2 Class

#imports
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)

from mos2class import mcell, safebandstructure, metric, xtohopvec
import tbplas as tb
import numpy as np
import math
import ast
import matplotlib.pyplot as plt
#-----------------------------------------------------


amax = 2
bmax = 2
cutoff_distance = 3
#------------------------------------------------------
#1.: Plotting e_min vs remaining hoppings
cell_1 = mcell("cell_1",0)
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
cell_2 = mcell("cell_2",0)
allhops_array, labels = cell_2.mhoppingtable()

#Comparison Cell 3:
cell_3 = mcell("cell_3",0)
cell_3_k_lenn, cell_3_bandss = cell_3.calcbands()


#Loop:
minhopa = 0.0
minhopb = 1.6
minhops = np.linspace(minhopa,minhopb,100)
metrics = []
for minhop in minhops:
    filtered_array = allhops_array[abs(allhops_array[:,6])>minhop]
    cell_2.changehops_toarr(filtered_array)
    metrics.append(metric(cell_2, cell_3_k_lenn,cell_3_bandss))
    print(metric(cell_2,cell_3_k_lenn,cell_3_bandss))

plt.gca().invert_xaxis()
plt.plot(nvec, metrics)
plt.ylim(top=3000) #Important!
plt.axhline(y=0,color='grey')
plt.xlabel("Number N of Hoppings / Total Number of Hoppings")
plt.ylabel("Error Metric m")
plt.title("Error Metric m for different N(E_min). \n amax = bmax = %s, cutoff-distance d = %s A" % (amax, cutoff_distance))
plt.savefig("pca_exp0_Metric_vsN_newmos2class.png")
#plt.show()
plt.clf()

#Save the data of the figure
mvsNfile = open("mvsN_fromlowEtohighE.txt",'w+')
mvsNfile.writelines(str(nvec)+"\n")
mvsNfile.writelines(str(metrics))
mvsNfile.close()



#------------------------------------------------------
#3.: Adding orbitals with E>Emin

#Experimental cell 3:
cell_4 = mcell("cell_4",0)
allhops_array, labels = cell_4.mhoppingtable()

#Comparison Cell 3:
cell_5 = mcell("cell_5",0)
cell_5_k_lenn, cell_5_bandss = cell_5.calcbands()


#Loop:
minhopa = 0.0
minhopb = 1.6
minhops = np.linspace(minhopa,minhopb,100)
metrics = []
for minhop in minhops:
    filtered_array = allhops_array[abs(allhops_array[:,6])>minhop]
    cell_4.changehops_toarr(filtered_array)
    metrics.append(metric(cell_4, cell_5_k_lenn, cell_5_bandss))
    print(metric(cell_4,cell_5_k_lenn, cell_5_bandss))

plt.gca().invert_xaxis()
plt.plot(minhops, metrics)
plt.xlabel("E_min")
plt.ylabel("Error Metric m")
plt.title("Error Metric m for different E_min.")
plt.savefig("pca_exp0_Metric_vsE_newmos2class.png")
#plt.show()
plt.clf()
