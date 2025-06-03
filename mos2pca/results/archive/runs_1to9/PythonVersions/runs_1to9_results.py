
import matplotlib.pyplot as plt
import tbplas as tb

#import mos2class
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from mos2class import mcell, safebandstructure, metric

import ast
#-----------------------------------------------------


ID = 9
E_min = 0.2 #Must be same as in pca_graddesc.py


#IDEAL CELL:
ideal_cell = mcell("Ideal",E_min)
idealhops = ideal_cell.mhoppings

#-----------------------------------------------------
#Functions:
#-----------------------------------------------------
def xtohopvec(x):
    hopvec = []
    for i in range(len(x)):
        if(idealhops[i]!=0):
            rn = idealhops[i][0]
            orb_i = idealhops[i][1]
            orb_j = idealhops[i][2]
            energy = idealhops[i][3]
            energy *= x[i]
            hopvec.append([rn,orb_i,orb_j,energy])
    return hopvec

def safe_Nvsm_graph(N_relative,m):
    name = "mvsN_fromlowEtohighE.txt"
    finalxfile = open(name,"r")
    content = []
    for line in finalxfile:
        content.append(line)
    nvec = ast.literal_eval(content[0])
    metrics = ast.literal_eval(content[1])

    #Plotting:
    plt.gca().invert_xaxis()
    plt.plot(nvec, metrics)
    plt.scatter(N_relative,m,color="Red")
    plt.xlabel("Number N of Hoppings / Total Number of Hoppings")
    plt.ylabel("Error Metric m")
    plt.title("Error Metric m for different hopping number N.")
    plt.savefig("graddesc_Nvsm_run"+str(ID)+".png")
    plt.clf()



#-----------------------------------------------------
#Main:
#-----------------------------------------------------

#Reading x from file
name = "graddesc_finalx_run" + str(ID) + ".txt"
finalxfile = open(name,"r")
content = finalxfile.read()
x = ast.literal_eval(content)

#Calculating hopvec:
hopvec = xtohopvec(x)
hopvec = [hop for hop in hopvec if abs(hop[3])>E_min]

#Build the final cell:
reducedhopcell = mcell("Mos2 Cell with grad-desc reduced hoppings", E_min)
reducedhopcell.changehops_tohopvec(hopvec)
N = reducedhopcell.mnhoppings

#Build the "true" cell:
truecell = mcell("Mos2 Cell with all hoppings", 0)
total = truecell.mnhoppings
m = metric(truecell,reducedhopcell)
safe_Nvsm_graph(N/total,m)

cellvector=[truecell,reducedhopcell]
filename="graddesc_bands_run"+str(ID)
title = "Bandstructures of MoS2 cell with all hoppings vs MoS2 cell with reduced hoppings"
safebandstructure(cellvector,filename,title)

indices = [i for i in range(len(x))]
plt.plot(indices,x,color = "Grey")
plt.scatter(indices,x,color="Red")
plt.xlabel("Index i")
plt.ylabel("Weight factor x[i]")
plt.title("Weight vector x")
plt.savefig("graddesc_xvector_run"+str(ID)+".png")
plt.clf()
