#Script for analyzing and plotting Results of Test runs of the Nesterov Gradient Descent PCA
#Based on the script "pca_graddesc.py"

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

path_runs = "finalruns/"
path_mvsN = "finalruns/mvsNanalysis/"

ID = 40 #ID of Run to evaluate
E_min = 0.1 #Must be same as in pca_graddesc.py


#IDEAL CELL:
ideal_cell = mcell("Ideal",E_min)
idealhops = ideal_cell.mhoppings
ideal_k_lenn, ideal_bandss = ideal_cell.calcbands()


#-----------------------------------------------------
#Functions:
#-----------------------------------------------------

def safe_Nvsm_graph(N_relative,m):
    #Get N vs M Information of Energyreduced Cell:
    name = path_mvsN + "mvsN_fromlowEtohighE.txt"
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
    plt.savefig(path_runs+"graddesc_Nvsm_run"+str(ID)+".png")
    plt.clf()



#-----------------------------------------------------
#Main:
#-----------------------------------------------------

#Reading x from file
name = path_runs + "graddesc_finalx_run" + str(ID) + ".txt"
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

#Plot Bandstructure Comparison with "true" cell
cellvector=[truecell,reducedhopcell]
filename= path_runs + "graddesc_bands_run"+str(ID)
title = "Bandstructures of MoS2 cell with all hoppings \n vs MoS2 cell with reduced hoppings, \n using Nesterov-Gradient-Descent"
safebandstructure(cellvector,filename,title)

#Plot Bandstructure Comparison with "true" and energyreduced cell
reducedenergycell = mcell("hoppings reduced in order of Energy",0)
allhops_array = reducedenergycell.mhoppingtable()[0]
allhops_array = allhops_array[allhops_array[:,6].argsort()[::-1]]
filtered_array = allhops_array[:N]
reducedenergycell.changehops_toarr(filtered_array)

cellvector=[truecell,reducedhopcell,reducedenergycell]
filename=path_runs + "graddesc_bands_redE_run"+str(ID)
title = "Bandstructures of MoS2 cell with all hoppings \n vs MoS2 cell with Nesterov-GD-reduced hoppings \n vs MoS2 cell with reduced hoppings in order of Energy"
safebandstructure(cellvector,filename,title)

#Plot x vector
indices = [i for i in range(len(x))]
plt.plot(indices,x,color = "Grey")
plt.scatter(indices,x,color="Red")
plt.xlabel("Index i")
plt.ylabel("Weight factor x[i]")
plt.title("Weight vector x")
plt.savefig(path_runs + "graddesc_xvector_run"+str(ID)+".png")
plt.clf()

#Plot hopping terms, reduced hopping Cell
filename=path_runs+"graddesc_cellplot_run"+str(ID)
reducedhopcell.mprimcell.plot(filename)

#Plot hopping terms, Comparison True Cell
#filename="graddesc_cellplot_comparison"
#truecell.mprimcell.plot(filename

#write down m and N_relative
mvsNfile = open(path_mvsN + "mvsN_run"+str(ID)+".txt", 'w+')
mvsNfile.writelines("["+str(N/total)+","+str(m)+"]")
mvsNfile.close()