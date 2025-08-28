#Plots m vs N graph for svd reductions

import numpy as np
import matplotlib.pyplot as plt
import ast

runpaths = []
runs = []
rundescriptions = []
runcolors = []
svdtypes = []
#-----------------------------------------------------
#PATHS, VARIABLES (IMPORTANT):
#-----------------------------------------------------

energyorderpath = "../resources/mvsN_energyorder/mvsN_energyorder.txt"
outputpath = "results/mvsn_svd" #without .png

#-----------------------------------------------------

setsize = 2

runpaths.append("results/newsvd/mvsNanalysis/")#runspath: path were mvsN_run___.txt files are stored
runs.append(range(1,23))
rundescriptions.append("new")
runcolors.append("green")
svdtypes.append("new")

runpaths.append("results/oldsvd/mvsNanalysis/")#runspath: path were mvsN_run___.txt files are stored
runs.append(range(1,5))
rundescriptions.append("old")
runcolors.append("blue")
svdtypes.append("old")



#-----------------------------------------------------
#MAIN:
#-----------------------------------------------------

def getnvsm(set):
    Nvec = []
    mvec = []
    Nvsm=[]
    for ID in runs[set]:
        stringid = svdtypes[set] + str(ID)
        name = runpaths[set] + "mvsN_run" + stringid + ".txt"
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
plt.plot(ns, metrics, color = "red", label="hoppings eliminated in ascending order of energy")

for i, rundescription in enumerate(rundescriptions):
    plt.plot(Nvsm[i][:,0], Nvsm[i][:,1], color = runcolors[i], label=rundescription)

#plt.ylim(top=500,bottom=0) #Important!
plt.xlabel("number N of hoppings / total number of hoppings")
plt.ylabel("error metric m in eV")
plt.legend(
    loc = 'upper center',
    bbox_to_anchor = (0.5,-0.18),
    ncol =1
)
plt.subplots_adjust(bottom=0.3)
#plt.title("error metric m for different number of hoppings")
plt.savefig(outputpath+".png")
plt.clf()