

import numpy as np
import tbplas as tb # Import the tbplas library
import matplotlib.pyplot as plt
from numpy.linalg import norm
from math import exp, sqrt
import io
import sys


#-----------------------------------------------------
#Functions:

#mhoppinglist(): returns a list of hopping informations for all hoppings in the testcell
def mhoppinglist(testcell):
    hoppinglist = []

    output = io.StringIO()
    sys.stdout = output
    testcell.print()
    sys.stdout = sys.__stdout__

    text = output.getvalue()
    text = text.split("terms:",1)[1]
    for line in text.splitlines():
        if(line == ""):
            continue
        for char in ("(",")",",","+","j"):
            line = line.replace(char, "%")
        linevec = line.split("%")
        linevec = [element for element in linevec if element.strip()]
        rn = [int(linevec[0]),int(linevec[1]),int(linevec[2])]
        orb_i = int(linevec[3])
        orb_j = int(linevec[4])
        energy = float(linevec[5])
        currenthop=[rn,orb_i,orb_j,energy]
        hoppinglist.append(currenthop)

    return hoppinglist


#-----------------------------------------------------
#CLASSES:
#mcell: represents a primitive cell with hoppings
#the custom classes and their variable names all start with an "m"

class mcell:
    
    def __init__(self,name,E_min):
        self.mname = name
        #Creating Primitive Cell
        self.mprimcell = tb.make_mos2_soc()
        self.mhoppings = mhoppinglist(self.mprimcell)
        self.mnhoppings = len(self.mhoppings)
        if (E_min > 1e-9):
            allhops_array = self.mhoppingtable()[0]
            filtered_array = allhops_array[abs(allhops_array[:,6])> E_min]
            self.changehops_toarr(filtered_array)

    #mclearhoppings deletes all hoppings of the mcell
    def mclearhoppings(self):
        for i in range(self.mnhoppings):
            self.mdelhop(0)
    
    #changehops_tohopvec clears all hoppings and adds all hoppings in newhoppingvector
    def changehops_tohopvec(self,newhoppingvector):
        self.mclearhoppings()
        self.mnewhoppings(newhoppingvector)

    #changehops_toarr clears all hoppings and adds hoppings based on numpy array
    def changehops_toarr(self,array):
        self.mclearhoppings()
        for k in range(array.shape[0]):
            rn = (array[k][1],array[k][2],array[k][3])
            i = array[k][4]
            j = array[k][5]
            e_hop = array[k][6]
            self.mnewhopping(rn,i,j,e_hop)

    #fundamental function for adding a hopping to the mcell class
    def mnewhopping(self, rn, i, j, e_hop):
        self.mprimcell.add_hopping(rn, i, j, e_hop)
        addedhop = [rn, i, j, e_hop]
        self.mhoppings.append(addedhop)
        self.mnhoppings = len(self.mhoppings)

    #adds all hoppings from vector "hoppings"
    def mnewhoppings(self, hoppings):
        for hop in hoppings:
            self.mprimcell.add_hopping(hop[0], hop[1], hop[2], hop[3]) #based on tbplas cell class
            self.mhoppings.append(hop)
        self.mnhoppings = len(self.mhoppings)

    #deletes hopping with index i
    def mdelhop(self,i):
        self.mprimcell.remove_hopping(self.mhoppings[i][0],self.mhoppings[i][1],self.mhoppings[i][2])
        del self.mhoppings[i]
        self.nhoppings = len(self.mhoppings)

    def mchangehop_energy(self, i, e_hop):
        self.mhoppings[i][3] = e_hop
        hop = self.mhoppings[i]
        self.mprimcell.remove_hopping(hop[0],hop[1],hop[2])
        self.mprimcell.add_hopping(hop[0], hop[1], hop[2], hop[3]) #based on tbplas cell class
    
    #returns the numpy array with the informations about all hoppings in the hoppingvector
    def mhoppingtable(self):
        vec = []
        header = "r,x,y,z,orb_i,orb_j,energy"
        for i in range(len(self.mhoppings)):
            currenthop = self.mhoppings[i]
            rn = currenthop[0]
            r = sqrt(rn[0]**2+rn[1]**2+rn[2]**2)
            currentdata = [r, rn[0], rn[1], rn[2], currenthop[1], currenthop[2], currenthop[3]]
            vec.append(currentdata)
        arr = np.array(vec)
        return arr, header


#Metric for the Difference between Bandstructures of cella and cellb
def metric(cella, cellb):
    cellvector = [cella.mprimcell,cellb.mprimcell]
    k_points = np.array([
        [0.0, 0.0, 0.0],    # Gamma
        [1./2, 0.0, 0.0],   # M
        [2./3, 1./3, 0.0],  # K
        [0.0, 0.0, 0.0],    # Gamma
    ])
    k_path, k_idx = tb.gen_kpath(k_points, [40, 40, 40])
    k_len=[]
    bands=[]
    numbands=[0,0]
    for i in range(2):
        k_lenn, bandss = cellvector[i].calc_bands(k_path,echo_details=False)
        k_len.append(k_lenn)
        bands.append(bandss)
        numbands[i] = bands[i].shape[1]

    if(numbands[0]!=numbands[1]):
        print("Metric Error: Different number of bands")
        return 0
    
    error = 0
    for i in range(numbands[0]):
        for j in range(bands[0].shape[0]):
            currenterror = abs(bands[0][j,i]-bands[1][j,i])
            if(currenterror >= 0.1):
                error += abs(bands[0][j,i]-bands[1][j,i])

    return error


#Safes Bandstructure Plot of all Cells in the mcellvector
def safebandstructure(mcellvector,filename,title):
    n= len(mcellvector)
    cellvector = [mcell.mprimcell for mcell in mcellvector]

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
                label = str(mcellvector[i].name) + "eV, " + str(mcellvector[i].mnhoppings) +" Hoppings"
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
    plt.title(title)
    pngname = filename + ".png"
    plt.savefig(pngname)
    #plt.show()
    plt.close()
