#prints orbitals in the R Roldán model to "orbitals.txt"

import tbplas as tb
from mos2class import mcell, mxtohopvec

cell = mcell("Test", 0)
cell.mclearhoppings_start()
orbitale = cell.mprimcell.orbitals
outputfile = open("orbitals.txt","w+")
for orbital in orbitale:
    outputfile.writelines(str(orbital)+"\n")
outputfile.close()