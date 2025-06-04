#Script for analyzing and plotting Results of Test runs of the Nesterov Gradient Descent PCA
#Based on the script "pca_graddesc.py"

import matplotlib.pyplot as plt
import tbplas as tb

#import mos2class
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from mos2class import mcell, msafebandstructure, mmetric, mxtohopvec

import ast
#-----------------------------------------------------
#VARIABLES (IMPORTANT):
pathsettings = "multilevelruns"

if (pathsettings == "finalruns"):
    path_runs = "finalruns/"
    path_output = "finalruns/formated_results/"
    path_fromlowEtohighEfile = "/home/bastian/bachelorarbeit_Projekte/Projekte/pca_mos2/ressources/mvsN_fromlowEtohighE.txt"
    path_mvsN_output = "finalruns/mvsNanalysis/"
elif (pathsettings == "multilevelruns"):
    path_runs = "multilevelruns/"
    path_output = "multilevelruns/formated_results/"
    path_fromlowEtohighEfile = "/home/bastian/bachelorarbeit_Projekte/Projekte/pca_mos2/ressources/mvsN_fromlowEtohighE.txt"
    path_mvsN_output = "multilevelruns/mvsNanalysis/"
elif (pathsettings == "customized"):
    path_runs = "highiterations/"
    path_output = "highiterations/formated_results/"
    path_fromlowEtohighEfile = "/home/bastian/bachelorarbeit_Projekte/Projekte/pca_mos2/ressources/mvsN_fromlowEtohighE.txt"
    path_mvsN_output = "highiterations/mvsNanalysis/"


ID = 122 #ID of Run to evaluate
E_min = 0.1 #Must be same as in pca_graddesc.py


#-----------------------------------------------------
#MAIN:

#IDEAL CELL:
ideal_cell = mcell("Ideal",E_min)
idealhops = ideal_cell.mhoppings
#ideal_bandss = ideal_cell.mcalcbands()


#-----------------------------------------------------
#Functions:
#-----------------------------------------------------

def safe_Nvsm_graph(N_relative,m):
    #Get N vs M Information of Energyreduced Cell:
    finalxfile = open(path_fromlowEtohighEfile,"r")
    content = []
    for line in finalxfile:
        content.append(line)
    nvec = ast.literal_eval(content[0])
    metrics = ast.literal_eval(content[1])

    #Plotting:
    plt.gca().invert_xaxis()
    plt.axhline(y=0,color='grey')
    plt.plot(nvec, metrics, color = "Red", label="Hoppings reduced in the order of increasing energies")
    plt.scatter(N_relative,m,color="Green", label="Hoppings reduced with Nesterov Gradient Descent")
    plt.ylim(top=1200,bottom=-300) #Important!
    plt.xlabel("Number N of Hoppings / Total Number of Hoppings")
    plt.ylabel("Error Metric m in eV")
    plt.title("Error Metric m for varying number N of hoppings. \n lambda_0 = ")
    plt.savefig(path_output+"graddesc_Nvsm_run"+str(ID)+".png")
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
hopvec = mxtohopvec(x,idealhops)
hopvec = [hop for hop in hopvec if abs(hop[3])>E_min]

#Build the final cell:
reducedhopcell = mcell("hoppings reduced by Nesterov grad-desc", E_min)
reducedhopcell.mchangehops_tohopvec(hopvec)
N = reducedhopcell.mnhoppings

#Build the "true" cell:
truecell = mcell("all hoppings", 0)
total = truecell.mnhoppings
true_bands = truecell.mcalcbands()

#Plot N vs M Graph
m = mmetric(reducedhopcell,true_bands)
safe_Nvsm_graph(N/total,m)

#Plot Bandstructure Comparison with "true" cell
cellvector=[truecell,reducedhopcell]
filename= path_output + "graddesc_bands_run"+str(ID)
title = "Bandstructures of MoS2 cell with all hoppings \n vs MoS2 cell with reduced hoppings, \n using Nesterov-Gradient-Descent"
msafebandstructure(cellvector,filename,title)

#Plot Bandstructure Comparison with "true" and energyreduced cell
reducedenergycell = mcell("hoppings reduced in order of Energy",0)
allhops_array = reducedenergycell.mhoppingtable()[0]
allhops_array = allhops_array[allhops_array[:,6].argsort()[::-1]]
filtered_array = allhops_array[:N]
reducedenergycell.mchangehops_toarr(filtered_array)

cellvector=[truecell,reducedhopcell,reducedenergycell]
filename=path_output + "graddesc_bands_redE_run"+str(ID)
title = "Bandstructures of MoS2 cell with all hoppings \n vs MoS2 cell with Nesterov-GD-reduced hoppings \n vs MoS2 cell with reduced hoppings in order of Energy"
msafebandstructure(cellvector,filename,title)

#Plot x vector
indices = [i for i in range(len(x))]
plt.plot(indices,x,color = "Grey")
plt.scatter(indices,x,color="Red")
plt.xlabel("Index i")
plt.ylabel("Weight factor x[i]")
plt.title("Weight vector x")
plt.savefig(path_output + "graddesc_xvector_run"+str(ID)+".png")
plt.clf()

#Plot hopping terms, reduced hopping Cell
filename=path_output+"graddesc_cellplot_run"+str(ID)
reducedhopcell.mprimcell.plot(filename)

#Plot hopping terms, Comparison True Cell
#filename="graddesc_cellplot_comparison"
#truecell.mprimcell.plot(filename

#write down m and N_relative
mvsNfile = open(path_mvsN_output + "mvsN_run"+str(ID)+".txt", 'w+')
mvsNfile.writelines("["+str(N/total)+","+str(m)+"]")
mvsNfile.close()

