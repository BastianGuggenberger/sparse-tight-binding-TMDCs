
import matplotlib.pyplot as plt
import tbplas as tb

#import mos2class
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from mos2class import mcell, safebandstructure, metric, xtohopvec

import ast
#-----------------------------------------------------


ID = 12
E_min = 0.2 #Must be same as in pca_graddesc.py


#IDEAL CELL:
ideal_cell = mcell("Ideal",E_min)
idealhops = ideal_cell.mhoppings
ideal_k_lenn, ideal_bandss = ideal_cell.calcbands()


#-----------------------------------------------------
#Functions:
#-----------------------------------------------------

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
    plt.axhline(y=0,color='grey')
    plt.plot(nvec, metrics)
    plt.scatter(N_relative,m,color="Red")
    plt.ylim(top=3000) #Important!
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
hopvec = xtohopvec(x,idealhops)
hopvec = [hop for hop in hopvec if abs(hop[3])>E_min]

#Build the final cell:
reducedhopcell = mcell("hoppings reduced by Nesterov grad-desc", E_min)
reducedhopcell.changehops_tohopvec(hopvec)
N = reducedhopcell.mnhoppings

#Build the "true" cell:
truecell = mcell("all hoppings", 0)
total = truecell.mnhoppings
true_k_lenn, true_bandss = truecell.calcbands()

#Plot N vs M Graph
m = metric(reducedhopcell,true_k_lenn,true_bandss)
safe_Nvsm_graph(N/total,m)

#Plot Bandstructure Comparison
cellvector=[truecell,reducedhopcell]
filename="graddesc_bands_run"+str(ID)
title = "Bandstructures of MoS2 cell with all hoppings \n vs MoS2 cell with reduced hoppings, \n using Nesterov-Gradient-Descent"
safebandstructure(cellvector,filename,title)

#Plot x vector
indices = [i for i in range(len(x))]
plt.plot(indices,x,color = "Grey")
plt.scatter(indices,x,color="Red")
plt.xlabel("Index i")
plt.ylabel("Weight factor x[i]")
plt.title("Weight vector x")
plt.savefig("graddesc_xvector_run"+str(ID)+".png")
plt.clf()

#Plot hopping terms, reduced hopping Cell
filename="graddesc_cellplot_run"+str(ID)
reducedhopcell.mprimcell.plot(filename)

#Plot hopping terms, Comparison True Cell
#filename="graddesc_cellplot_comparison"
#truecell.mprimcell.plot(filename