import matplotlib.pyplot as plt
import tbplas as tb
from mos2class import mcell, msafebandstructure, mmetric, mxtohopvec

truecell = mcell("true",0)
total=truecell.mnhoppings

truebands = truecell.mcalcbands()

eminvec = [0.01 * i for i in range(50)]
for emin in eminvec:
    testcell = mcell("test",emin)
    output = "totalhops: " + str(testcell.mnhoppings) + ",  relhops: " + str(testcell.mnhoppings/total) + ",  E_min: " + str(emin) + ",  metric: " + str(mmetric(testcell,truebands))
    print(output)
