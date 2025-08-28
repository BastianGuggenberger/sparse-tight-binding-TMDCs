#Script for analyzing and plotting Results of SVD based hopping reductions

import matplotlib.pyplot as plt
import tbplas as tb
from mos2class import mcell, msafebandstructure, mmetric, mxtohopvec

import ast
#-----------------------------------------------------
#VARIABLES (IMPORTANT):
IDset = range(1,23) #IDs of Runs to evaluate
E_min = 0.1 #Must be same as in pca_graddesc.py
svdtype = "new"

if(svdtype == "old"):
    path_runs = "results/oldsvd/"
    path_output = "results/oldsvd/formated_results/"
    path_mvsN_output = "results/oldsvd/mvsNanalysis/"
elif(svdtype == "new"):
    path_runs = "results/newsvd/"
    path_output = "results/newsvd/formated_results/"
    path_mvsN_output = "results/newsvd/mvsNanalysis/"
    
path_fromlowEtohighEfile = "../resources/mvsN_energyorder/mvsN_energyorder.txt"


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
    fromlowetohighEfile = open(path_fromlowEtohighEfile,"r")
    content = []
    for line in fromlowetohighEfile:
        content.append(line)
    nvec = ast.literal_eval(content[0])
    metrics = ast.literal_eval(content[1])

    #Plotting:
    plt.gca().invert_xaxis()
    plt.axhline(y=0,color='grey')
    plt.plot(nvec, metrics, color = "Red", label="Hoppings reduced in the order of increasing energies")
    plt.scatter(N_relative,m,color="Green", label="Hoppings reduced with Nesterov Gradient Descent")
    plt.ylim(top=500,bottom=0) #Important!
    plt.xlabel("Number N of Hoppings / Total Number of Hoppings")
    plt.ylabel("Error Metric m in eV")
    plt.title("Error Metric m for varying number N of hoppings. \n lambda_0 = ")
    plt.savefig(path_output+"graddesc_Nvsm_run"+idstring+".png")
    plt.clf()



#-----------------------------------------------------
#Main:
#-----------------------------------------------------

for ID in IDset:
    idstring = svdtype + str(ID)

    #Reading hopvec from file
    name = path_runs + "graddesc_finalhopvec_" + idstring + ".txt"
    finalhopsfile = open(name,"r")
    content = finalhopsfile.read()
    hopvec = ast.literal_eval(content)

    #Build the final cell:
    reducedhopcell = mcell("hoppings reduced by svd", E_min)
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
    filename= path_output + "graddesc_bands_run"+ idstring
    title = "Bandstructures of MoS2 cell with all hoppings \n vs MoS2 cell with reduced hoppings, \n using Nesterov-Gradient-Descent"
    msafebandstructure(cellvector,filename,title)

    #Plot Bandstructure Comparison with "true" and energyreduced cell
    reducedenergycell = mcell("hoppings reduced in order of Energy",0)
    allhops_array = reducedenergycell.mhoppingtable()[0]
    allhops_array = allhops_array[allhops_array[:,6].argsort()[::-1]]
    filtered_array = allhops_array[:N]
    reducedenergycell.mchangehops_toarr(filtered_array)

    cellvector=[truecell,reducedhopcell,reducedenergycell]
    filename=path_output + "graddesc_bands_redE_run"+ idstring
    title = "Bandstructures of MoS2 cell with all hoppings \n vs MoS2 cell with Nesterov-GD-reduced hoppings \n vs MoS2 cell with reduced hoppings in order of Energy"
    msafebandstructure(cellvector,filename,title)


    #Plot hopping terms, reduced hopping Cell
    filename=path_output+"graddesc_cellplot_run"+str(ID)
    reducedhopcell.mprimcell.plot(filename)

    #Plot hopping terms, Comparison True Cell
    #filename="graddesc_cellplot_comparison"
    #truecell.mprimcell.plot(filename

    #write down m and N_relative
    mvsNfile = open(path_mvsN_output + "mvsN_run"+idstring+".txt", 'w+')
    mvsNfile.writelines("["+str(N/total)+","+str(m)+"]")
    mvsNfile.close()
    print("Safed results of run " + idstring)
print("Safed all results succesfully")
