#Script for analysing the N vs m history of the gradient descent algorithm

import matplotlib.pyplot as plt
import numpy as np

#-----------------------------------------------------
#VARIABLES (IMPORTANT):
path_runs = "../../highiterations_2.0/"
path_output = "../../importantresults/"
ID = 1007 #ID of Run to evaluate

#-----------------------------------------------------

name = path_runs + "graddesc_history_run" + str(ID) + ".txt"
historyfile = open(name,"r")

Nvector_seperated = [[],[],[]] #3 vectors for orders 0,1,2
Nvector = []
mvector = []

for line in historyfile:
    if(line[0]=="N"):
        if(line[6]=="1"):
            linesep = line.split("=")
            content = linesep[1].split("/")
            Nvector_seperated[1].append(int(content[0]))
            continue
        elif(line[6]=="2"):
            linesep = line.split("=")
            content = linesep[1].split("/")
            Nvector_seperated[2].append(int(content[0]))
        else:
            print("Error: could not recognize order.")

    if(line[0]=="m"):
        linesep = line.split("=")
        mvector.append(float(linesep[1]))

Nvector = [(N_1 + N_2) for N_1, N_2 in zip(Nvector_seperated[1],Nvector_seperated[2])]

if(len(Nvector)!= len(mvector)):
    print("Error: N and m informations do not correspond in length.")

indices = [i for i in range(len(Nvector))]
ratio = []
for i in range(len(Nvector)):
    if(mvector[i]==0):
        ratio.append(0)
    else:
        ratio.append(40*float(Nvector[i])/mvector[i])


#Plot alllevel N vs m history
plt.plot(indices,Nvector, color ="violet", label="N(x[i]~0), total")
plt.plot(indices,mvector, color ="black", label = "metric m")
#plt.plot(indices,ratio,color="green", label = "m/N ratio")
plt.xlabel("NGD iterations")
plt.ylabel("N(x[i]~0), metric m in eV")
#plt.ylim(top = 400, bottom = 0)
plt.legend()

name_Nalllevels = path_output + "historyanalysis_alllevels_run" +str(ID)+ ".png"
plt.savefig(name_Nalllevels)
plt.clf()


#Plot seperated level N vs m history
plt.plot(indices,mvector, color ="black", label = "metric m",linewidth = 0.2)
plt.plot(indices,Nvector, color ="violet", label="N(x[i]~0), total",linewidth = 0.8)
plt.plot(indices,Nvector_seperated[1], color ="red", label="N(x[i]~0), order 1",linewidth = 1.2)
plt.plot(indices,Nvector_seperated[2], color ="blue", label="N(x[i]~0), order 2",linewidth = 1.2)
#plt.plot(indices,ratio,color="green", label = "m/N ratio")
#plt.ylim(top = 400, bottom = 0)
plt.xlabel("NGD iterations")
plt.ylabel("N(x[i]~0), metric m in eV")
plt.legend()

name_Nseplevels = path_output + "historyanalysis_seperatedlevels_run" +str(ID)+ ".png"
plt.savefig(name_Nseplevels)
