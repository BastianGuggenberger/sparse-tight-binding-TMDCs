#Compares the bandstructures of SVD based reduced MoS2 cell
#iterates for different n_components

#import class mcell and imports from mos2class:
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from mos2class  import mcell, msafebandstructure

#-----------------------------------------------------------------------
E_min = 0.1

for N in range(2,5):
    idealcell = mcell("Ideal Cell",E_min)
    svdcell = mcell("SVD based hoppings reduced Cell",E_min)
    print("Original Hoppings: ", svdcell.mnhoppings)

    svdcell.msvdtruncate(N)
    print("New Hoppings: ",svdcell.mnhoppings)
    path = "finalhopvec_"+str(N)+".txt"
    finalhopfile = open(path, "w+")
    finalhopfile.write(str(svdcell.mhoppings.copy()))
    finalhopfile.close()

    cellvec = [idealcell,svdcell]
    msafebandstructure(cellvec,"comparison_N=" + str(N), "Comparison of SVD reduced vs Ideal Hoppings, N = " + str(N))