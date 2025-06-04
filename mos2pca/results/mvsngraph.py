#Plots m vs N graph

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
energyorderpath = "/home/bastian/bachelorarbeit_Projekte/Projekte/pca_mos2/ressources/mvsN_fromlowEtohighE.txt"
outputpath = "multilevelruns/mvsNanalysis/mvsNcomparison_102_112_andallhighiterationsruns_andallfinalruns" #without .png

#-----------------------------------------------------
#sets:
setsize = 3 #how many different paths (sets) are added

#all paths, runs and description of the sets we want to analyze
runpaths.append("multilevelruns/mvsNanalysis/")#runspath: path were mvsN_run___.txt files are stored
runs.append([102,112])
rundescriptions.append("Hoppings reduced with Multilevel-NGD, 1200 Iterations")
runcolors.append("green")

runpaths.append("highiterations/mvsNanalysis/")#runspath: path were mvsN_run___.txt files are stored
runs.append([54,60,61,62,63])
rundescriptions.append("Hoppings reduced with NGD, 1200 Iterations")
runcolors.append("blue")

runpaths.append("finalruns/mvsNanalysis/")#runspath: path were mvsN_run___.txt files are stored
runs.append(range(40,54))
rundescriptions.append("Hoppings reduced with NGD, 400 Iterations")
runcolors.append("orange")


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
    content.append(line)
energyorderfile.close()
ns = ast.literal_eval(content[0])
metrics = ast.literal_eval(content[1])

#Plotting:
plt.gca().invert_xaxis()
plt.axhline(y=0,color='grey')
plt.plot(ns, metrics, color = "red", label="Hoppings reduced in the order of increasing energies")

for i, rundescription in enumerate(rundescriptions):
    plt.plot(Nvsm[i][:,0], Nvsm[i][:,1], color = runcolors[i], label=rundescription)

plt.ylim(top=1200,bottom=-500) #Important!
plt.xlabel("Number N of Hoppings / Total Number of Hoppings")
plt.ylabel("Error Metric m in eV")
plt.legend()
plt.title("Error Metric m for varying number N of hoppings.")
plt.savefig(outputpath+".png")
plt.clf()