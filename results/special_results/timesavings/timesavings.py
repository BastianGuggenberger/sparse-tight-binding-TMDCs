#Provides the timesavings vs error plot

import time
import matplotlib.pyplot as plt
import tbplas as tb
from mos2class import mcell, msafebandstructure, mmetric, mxtohopvec, mmaxerroratbandgap, mtotalerrorperkpoint

import ast
#-----------------------------------------------------
#VARIABLES (IMPORTANT):
#-----------------------------------------------------
path_runs = "../../highiterations_newk/"
IDs = range(2100,2133)
N_IDs = len(IDs)
E_min = 0.1


#-----------------------------------------------------
#MAIN:
#-----------------------------------------------------
#IDEAL CELL:
ideal_cell = mcell("Ideal",E_min)
idealhops = ideal_cell.mhoppings.copy()


originaltime = 0
timesvec = [0 for i in range(N_IDs)]
errorsvec = [0 for i in range(N_IDs)]

for i in range(10):

    start = time.time()
    ideal_bands = ideal_cell.mcalcbands()
    end = time.time()
    originaltime += end-start

    for j, ID in enumerate(IDs):
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
        start = time.time()
        reducedhopcell.mcalcbands()
        end = time.time()

        timesvec[j]+=(end-start)
        errorsvec[j]+=(mtotalerrorperkpoint(reducedhopcell,ideal_bands))

#calculate means:

originaltime=originaltime/10
timesvec = [totaltime/10 for totaltime in timesvec]
errorsvec = [totalerror/10 for totalerror in errorsvec]

timesvec.append(originaltime)
errorsvec.append(0.0)

timesvec = [calctime/originaltime for calctime in timesvec]

plt.scatter(timesvec,errorsvec)
plt.xlabel("reduced BS calculation time / original BS calculation time")
plt.ylabel("accuracy loss: bandstructure error per band and k-point in eV")
plt.savefig("timesavings",dpi=300)