#Class for adding, deleting, analyzing mos2 hopping terms with different methods

#Based on tbplas method make_mos2_soc, which is an implementation of:
    # R Roldán et al 2014 2D Mater. 1 034003
    # https://www.tbplas.net/_api/tbplas.make_mos2_soc.html

#the custom class variable-, method- and all custom function- names start with an "m" for better distinguishability to tbplas


import numpy as np
import tbplas as tb # Import the tbplas library
from math import exp, sqrt
import io
import sys
import ast
from scipy.linalg import svd
import time
#-----------------------------------------------------
#PATHS:
idealhoppath= "/home/bastian/bachelorarbeit_Projekte/Projekte/pca_mos2/ressources/idealhops.txt"

#-----------------------------------------------------
#Global Variables:

k_points = np.array([
        [0.0, 0.0, 0.0],    # Gamma
        [1./2, 0.0, 0.0],   # M
        [2./3, 1./3, 0.0],  # K
        [0.0, 0.0, 0.0],    # Gamma
    ])

k_path, k_idx = tb.gen_kpath(k_points, [40, 40, 40])
k_path_efficient, k_idx_efficient = tb.gen_kpath(k_points, [15, 15, 15]) #kpath with less points for faster calculation
k_label = ["G", "M", "K", "G"]

colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628"]

#idealhoppingslist: list of all hoppings in R Roldán et al paper.
#list is built with hoplistbuiler.py
#imaginary part is not taken into account
idealhoppingslist = []
idealhopfile = open(idealhoppath,"r")
for line in idealhopfile:
    vec = ast.literal_eval(line)
    idealhoppingslist.append([vec[0],int(vec[1]),int(vec[2]),vec[3]])


#-----------------------------------------------------
#CLASSES:

#mcell: represents a primitive cell with hoppings
class mcell:
    
    def __init__(self,name,E_min):
        self.mname = name
        self.mnorb = 22 #max orb is 21 - there are 22 different orbitals in the R Roldán et al paper.
        self.mE_min = E_min #Minimum Energy of a hopping to be evaluated

        #Creating Primitive Cell
        self.mprimcell = tb.make_mos2_soc()

        #Changing hoppings of primitive cell to hoppings in idealhoppinglist - relevant, because make_mos2_soc() double-counts some hoppings and has imaginary parts (soc)
        self.mhoppings = idealhoppingslist.copy()
        self.mnhoppings = len(self.mhoppings)
        self.mchangehops_tohopvec(idealhoppingslist.copy())
        
        #Filtering out hoppings with |E| < Emin
        allhops_array = self.mhoppingtable()[0]
        filtered_array = allhops_array[abs(allhops_array[:,6])> E_min]
        self.mchangehops_toarr(filtered_array)

    #-----------------------------------------------------
    #METHODS of class mcell:

    #deletes all hoppings of the mcell
    def mclearhoppings(self):
        for i in range(self.mnhoppings):
            self.mdelhop(0)
    
    #clears all hoppings and adds all hoppings in newhoppingvector
    def mchangehops_tohopvec(self,newhoppingvector):
        self.mclearhoppings()
        self.mnewhoppings(newhoppingvector)

    #clears all hoppings and adds hoppings based on numpy array
    def mchangehops_toarr(self,array):
        self.mclearhoppings()
        for k in range(array.shape[0]):
            rn = (array[k][1],array[k][2],array[k][3])
            i = int(array[k][4])
            j = int(array[k][5])
            e_hop = array[k][6]
            self.mnewhopping(rn,i,j,e_hop)

    #fundamental method for adding a hopping to the mcell class
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
        self.mnhoppings = len(self.mhoppings)

    #changes energy of ith hopping to e_hop
    def mchangehop_energy(self, i, e_hop):
        self.mhoppings[i][3] = e_hop
        hop = self.mhoppings[i].copy()
        self.mprimcell.add_hopping(hop[0], hop[1], hop[2], hop[3]) #tbplas add_hopping can also update an existing hop

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

    #calculates the bandstructure of the cell, "efficient" calculates less points, for faster calculation
    def mcalcbands(self, efficient = False):
        if (efficient):
            path = k_path_efficient
        else:
            path = k_path
        klen, bands = self.mprimcell.calc_bands(path,echo_details=False)
        return bands
    
    #method for getting the tbplas internal hopping index of a given hop
    #the index of tbplas hopping does not usually correspond with the mcell.mhoppings index
    def mget_hopping_index(self,hop):
        hop_ind = np.array([hop[0][0],hop[0][1],hop[0][2],hop[1],hop[2]])
        index = None
        for i, current_hop_ind in enumerate(self.mprimcell.hop_ind):
            if(np.array_equal(hop_ind,current_hop_ind)):
                index = i

        return index
    
    #analyzes the neighbours degrees of all hoppings in the cell
    def mget_neighbours_orders(self):
        hopvec = self.mhoppings.copy()
        dset = [0.0, 0.24, 0.32]
        dvec = []
        for hop in hopvec:
            index = self.mget_hopping_index(hop)
            r = self.mprimcell.dr_nm[index]
            d = np.linalg.norm(r)
            d = np.round(d,2)
            dvec.append(d)
            if d not in dset:
                print("unknown hopping distance d=", d)
                return None

        #Categorize
        n = [0,0,0]
        for d in dvec:
            uncategorized = True
            for i, compd in enumerate(dset):
                if(d==compd):
                    uncategorized = False
                    n[i]+=1
            if(uncategorized):
                print("Hop not categorized! d=", d)
                return None
            
        return n
    
    #reduces hoppings based on SVD and truncation
    def msvdtruncate(self,n_components):
        #1. Get T - Matrix
        hopvec = self.mhoppings.copy()
        hopdict = {}
        for hop in hopvec:
            if(hop[0] in hopdict.keys()):
                hopdict[hop[0]][hop[1],hop[2]] = hop[3]
                hopdict[hop[0]][hop[2],hop[1]] = hop[3].conj()

            else:
                hopdict[hop[0]] = np.zeros((self.mnorb,self.mnorb))
                hopdict[hop[0]][hop[1],hop[2]] = hop[3]
                hopdict[hop[0]][hop[2],hop[1]] = hop[3].conj()

        r_vec = list(hopdict.keys())
        hop_vec = list(hopdict.values())
        
        T = np.array([hop.flatten() for hop in hop_vec])

        #2. Perform SVD
        U, s, Vh = svd(T, full_matrices=False)

        #3. Truncate
        U = U[:, :n_components]
        s = s[:n_components]
        Vh = Vh[:n_components, :]
        T = U @ np.diag(s) @ Vh
                    
        #4. Calculate New hopping list
        newhopvec = []
        for r_index,R in enumerate(r_vec):
            hamiltonian = T[r_index].reshape(self.mnorb,self.mnorb)
           
            #Enforce Hermicity of Hamiltonian:#(Clears Imaginary part, only compatible with no SOC)
            hamiltonian = 0.5*(hamiltonian + hamiltonian.conj().T)

            for i in range(self.mnorb):
                for j in range(i+1):
                    Ehop = hamiltonian[i,j]
                    if(abs(Ehop)>self.mE_min):
                        newhopvec.append([R,i,j,Ehop])   
        
        self.mchangehops_tohopvec(newhopvec)  




#-----------------------------------------------------
#FUNCTIONS:

#Metric for the Difference between Bandstructures of cella and cellb, for cellb the bands have to be given
#"efficient" only compares less points for faster calculation
def mmetric(cell, comparison_bands, efficient = False):

    #"efficient": only 45 instead of 120 k-Points are calculated
    if (efficient):
        path = k_path_efficient
    else:
        path = k_path

    k_len, bands = cell.mprimcell.calc_bands(path,echo_details=False) #line takes about 80% of the total time for an iteration in pca_graddesc.py

    bands_vector=[bands,comparison_bands]
    numbands_vector = [bands_vector[0].shape[1],bands_vector[1].shape[1]]

    if(bands_vector[0].shape[0]!=bands_vector[1].shape[0]):
        print("Metric Error: Different number of k-Points")
        return 0
    if(numbands_vector[0]!=numbands_vector[1]):
        print("Metric Error: Different number of bands")
        return 0
    

    #Error calculation:
    error = 0
    for i in range(numbands_vector[0]):
        for j in range(bands_vector[0].shape[0]):
            currenterror = (bands_vector[0][j,i]-bands_vector[1][j,i])**2
            if(currenterror >= 0.1):
                error += currenterror

    return error


#Safes Bandstructure Plot of all Cells in the mcellvector
def msafebandstructure(mcellvector,filename,title):
    import matplotlib.pyplot as plt

    n= len(mcellvector)
    cellvector = [mcell.mprimcell for mcell in mcellvector]

    k_len_vector=[]
    bands_vector=[]
    for i in range(n):
        k_len, bands = cellvector[i].calc_bands(k_path)
        k_len_vector.append(k_len)
        bands_vector.append(bands)

    #Plotting the Bandstructure
    for i in range(n):
        color = colors[i]
        num_bands = bands_vector[i].shape[1]
        width = 0.7
        for j in range(num_bands):
            if(j==0):
                label = str(mcellvector[i].mname) + " , " + str(mcellvector[i].mnhoppings) +" Hoppings"
                plt.plot(k_len_vector[i], bands_vector[i][:, j], color=color, linewidth=width, label = label)
            else:
                plt.plot(k_len_vector[i], bands_vector[i][:, j], color=color, linewidth=width)

        for idx in k_idx:
            plt.axvline(k_len_vector[i][idx], color='k', linewidth=1.0)
        plt.xlim((0, np.amax(k_len_vector[i])))
        plt.xticks(k_len_vector[i][k_idx], k_label)
    plt.xlabel("k (1/nm)")
    plt.ylabel("Energy (eV)")
    #plt.tight_layout()
    plt.legend()
    plt.title(title)
    pngname = filename + ".png"
    plt.tight_layout()
    plt.savefig(pngname)
    #plt.show()
    plt.close()


#Get the Hopping vector from the weight vector x
def mxtohopvec(x,idealhops):
    hopvec = []
    for i in range(len(x)):
        if(idealhops[i]!=0):
            rn = idealhops[i][0]
            orb_i = idealhops[i][1]
            orb_j = idealhops[i][2]
            energy = idealhops[i][3]
            energy *= x[i]
            hopvec.append([rn,orb_i,orb_j,energy])
    return hopvec