#converts idealhops.txt to multilevelhops.txt
#idealhops.txt is not sorted into orders of neighbourhood
#multilevelhops.txt is sorted into oders of neighbourhood

import numpy as np
from mos2class.mos2class import mcell, msafebandstructure, mmetric, mxtohopvec

#-----------------------------------------------------
#PATHS:
idealhoppath= "../idealhoplist/idealhoplist.txt"
outputpath = "idealhoplist_multilevel.txt"

#-----------------------------------------------------


#The following code largely corresponds to the method "mgetneighboursorders" of class mcell.

#get all roldan hoppings
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
print("Succesfully stored hoppings in 'idealhoplist_multilevel.txt'")
