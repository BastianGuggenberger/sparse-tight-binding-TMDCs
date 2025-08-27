import tbplas as tb

#import mos2class
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from Projekte.reduced_mos2.resources.mos2class.mos2class import mcell, mxtohopvec

cell = mcell("Test", 0)
cell.mclearhoppings_start()
orbitale = cell.mprimcell.orbitals
outputfile = open("orbitals.txt","w+")
for orbital in orbitale:
    outputfile.writelines(str(orbital)+"\n")
outputfile.close()