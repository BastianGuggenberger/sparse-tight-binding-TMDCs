#Script for analysing the N vs m history of the gradient descent algorithm

import matplotlib.pyplot as plt
import numpy as np

#-----------------------------------------------------
#VARIABLES (IMPORTANT):
path_runs = "finalruns/"
path_output = "finalruns/formated_results/"
ID = 40 #ID of Run to evaluate

#-----------------------------------------------------

name = path_runs + "graddesc_history_run" + str(ID) + ".txt"
historyfile = open(name,"r")

Nvector = []
mvector = []

for line in historyfile:
    if(line[0]=="n"):
        continue
    if(line[0]=="N"):
        linesep = line.split("=")
        Nvector.append(int(linesep[1]))
        continue
    if(line[0]=="m"):
        linesep = line.split("=")
        mvector.append(float(linesep[1]))

if(len(Nvector)!= len(mvector)):
    print("Error: N and m informations do not correspond in length.")
indices = [i for i in range(len(Nvector))]
ratio = []
for i in range(len(Nvector)):
    if(mvector[i]==0):
        ratio.append(0)
    else:
        ratio.append(40*float(Nvector[i])/mvector[i])

plt.plot(indices,Nvector, color ="violet",label = "N(x[i]~0)")
plt.plot(indices,mvector, color ="black",label = "metric m")
plt.plot(indices,ratio,color="green", label = "m/N ratio")
plt.xlabel("gd iterations")
plt.ylabel("N(x[i]~0), metric m in eV, N/m ratio (scaled)")
plt.legend()

name_Nalllevels = path_output + "historyanalysis_alllevels_run" + str(ID) +".png"
plt.savefig(name_Nalllevels)