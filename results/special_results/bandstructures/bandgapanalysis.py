#Prints the maximum error of the reduced model near the bandgap.

import matplotlib.pyplot as plt
import tbplas as tb
from mos2class import mcell, msafebandstructure, mmetric, mxtohopvec, mmaxerroratbandgap

import ast
#-----------------------------------------------------
#VARIABLES (IMPORTANT):
#-----------------------------------------------------
path_runs = "../../bandgapruns_newk/"
ID = 2413
E_min = 0.1

#-----------------------------------------------------
#MAIN:
#-----------------------------------------------------
#IDEAL CELL:
ideal_cell = mcell("Ideal",E_min)
idealhops = ideal_cell.mhoppings.copy()
ideal_bands = ideal_cell.mcalcbands()


#Reading x from file
name = path_runs + "graddesc_finalx_run" + str(ID) + ".txt"
finalxfile = open(name,"r")
content = finalxfile.read()
content = content.replace("np.float64", "")
x = ast.literal_eval(content)

#Calculating hopvec:
hopvec = mxtohopvec(x,idealhops)
hopvec = [hop for hop in hopvec if abs(hop[3])>E_min]

#Build the final cell:
reducedhopcell = mcell("hoppings reduced by Nesterov grad-desc", E_min)
reducedhopcell.mchangehops_tohopvec(hopvec)
N = reducedhopcell.mnhoppings

maxerroratbandgap = mmaxerroratbandgap(reducedhopcell,ideal_bands)
print(maxerroratbandgap)