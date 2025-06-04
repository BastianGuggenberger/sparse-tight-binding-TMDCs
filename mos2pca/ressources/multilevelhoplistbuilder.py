#converts idealhops.txt to multilevelhops.txt
#idealhops.txt is not sorted into orders of neighbourhood
#multilevelhops.txt is sorted into oders of neighbourhood

import numpy as np
import ast

#import mos2class:
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from mos2class import mcell, msafebandstructure, mmetric, mxtohopvec

#-----------------------------------------------------
#PATHS:
idealhoppath= "/home/bastian/bachelorarbeit_Projekte/Projekte/pca_mos2/ressources/idealhops.txt"
outputpath = "/home/bastian/bachelorarbeit_Projekte/Projekte/pca_mos2/ressources/multilevelhops.txt"

#-----------------------------------------------------


#The following code largely corresponds to the method "mgetneighboursorders" of class mcell.

idealcell = mcell("Ideal Cell", 0)

hopvec = idealcell.mhoppings.copy()

dset = [0.0, 0.24, 0.32]
ml_hopvec = [[],[],[]]

for hop in hopvec:
    index = idealcell.mget_hopping_index(hop)
    r = idealcell.mprimcell.dr_nm[index]
    d = np.linalg.norm(r)
    d = np.round(d,2)

    #Categorize hopping:
    if d not in dset:
        print("unknown hopping distance d=", d)
    
    uncategorized = True
    for i, compd in enumerate(dset):
        if(d==compd):
            uncategorized = False
            ml_hopvec[i].append(hop)
    if(uncategorized):
        print("Hop not categorized! d=", d)

outputfile = open(outputpath,"w+")
outputfile.write(str(ml_hopvec))
print("Succesfully stored hoppings in 'multilevelhops.txt'")
