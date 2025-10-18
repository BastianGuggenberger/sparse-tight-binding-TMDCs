#The scripts take the mvsN (metric vs remaining number of hoppings) information of each run from "reduced_mos2/results/-runtype-/mvsNanalysis/mvsN_run1234.txt" and stores a plot of this information in the folder "reduced_mos2/results/-runtype-/formated_results/".  
#It also plots a comparison of all mvsNs of all different methods and stores this in the "importantresults" folder.
import numpy as np
import matplotlib.pyplot as plt
import ast

runpaths = []
runs = []
rundescriptions = []
runcolors = []
#-----------------------------------------------------
#PATHS, VARIABLES (IMPORTANT):
#-----------------------------------------------------

#inputpath: path were mvsN_fromlowEtohighE.txt is stored
energyorderpath = "../../../resources/mvsN_energyorder/mvsN_energyorder.txt"
outputpath = "../../importantresults/mvsNcomparison_mlvsal" #without .png

#-----------------------------------------------------
#sets:
setsize = 2 #how many different paths (sets) are added
"""
#all paths, runs and description of the sets we want to analyze
runpaths.append("../../multilevelruns/mvsNanalysis/")#runspath: path were mvsN_run___.txt files are stored
multilevelruns = [122,123,125,126]
for moreruns in range(127,156):
    multilevelruns.append(moreruns)
runs.append(multilevelruns)
rundescriptions.append("hoppings reduced with Multilevel-NGD, 1200 Iterations")
runcolors.append("green")
"""

"""
runpaths.append("../../multilevelruns_2.0/mvsNanalysis/")#runspath: path were mvsN_run___.txt files are stored
runs.append(range(1200,1232))
rundescriptions.append("hoppings reduced with separated-orders-NGD, 1200 Iterations")
runcolors.append("green")
"""


runpaths.append("../../lowiterations_newk/mvsNanalysis/")#runspath: path were mvsN_run___.txt files are stored
runs.append(range(2000,2033))
rundescriptions.append("hoppings reduced with NGD, 400 iterations")
runcolors.append("grey") 

runpaths.append("../../highiterations_newk/mvsNanalysis/")#runspath: path were mvsN_run___.txt files are stored
runs.append(range(2100,2133))
rundescriptions.append("hoppings reduced with NGD, 1200 iterations")
runcolors.append("blue")

"""
runpaths.append("../../finalruns/mvsNanalysis/")#runspath: path were mvsN_run___.txt files are stored
runs.append(range(40,54))
rundescriptions.append("hoppings reduced with NGD, 400 Iterations")
runcolors.append("orange")
"""



"""
runpaths.append("../../startzeroruns/mvsNanalysis/")#runspath: path were mvsN_run___.txt files are stored
startzeroruns = range(400,409)
runs.append(startzeroruns)
rundescriptions.append("Hoppings reduced with NGD, 400 Iterations, x_0 = 0")
runcolors.append("orange")
"""

#-----------------------------------------------------
#MAIN:
#-----------------------------------------------------

def getnvsm(set):
    Nvec = []
    mvec = []
    Nvsm=[]
    for ID in runs[set]:
        name = runpaths[set] + "mvsN_run" + str(ID) + ".txt"
        finalxfile = open(name,"r")
        content = finalxfile.read()
        x = ast.literal_eval(content)
        finalxfile.close()
        Nvsm.append(x)
        Nvec.append(x[0])
        mvec.append(x[1])

    #Sorting
    Nvsm = np.array(Nvsm)
    Nvsm = Nvsm[Nvsm[:,0].argsort()]
    return Nvsm

Nvsm = [None] *setsize
for i in range(setsize):
    Nvsm[i] = getnvsm(i)


#Get N vs M Information of Energyreduced Cell:
energyorderfile = open(energyorderpath,"r")
content = []
for line in energyorderfile:
    line = line.replace("np.float64", "")
    content.append(line)
energyorderfile.close()
ns = ast.literal_eval(content[0])
metrics = ast.literal_eval(content[1])

#Plotting:
plt.gca().invert_xaxis()
plt.axhline(y=0,color='grey')
plt.plot(ns, metrics, color = "red", label="hoppings eliminated in ascending order of energy")

for i, rundescription in enumerate(rundescriptions):
    plt.plot(Nvsm[i][:,0], Nvsm[i][:,1], color = runcolors[i], label=rundescription)

plt.ylim(top=500,bottom=0) #Important!
plt.xlabel("number N of hoppings / total number of hoppings")
plt.ylabel("error metric m in (eV)^2")
plt.legend(
    loc = 'upper center',
    bbox_to_anchor = (0.5,-0.18),
    ncol =1
)
plt.subplots_adjust(bottom=0.3)
#plt.title("error metric m for different number of hoppings")
plt.savefig(outputpath+".png")
plt.clf()