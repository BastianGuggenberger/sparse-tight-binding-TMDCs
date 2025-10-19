#Creates a plot of a bandstructure comparison for evaluating the performance of the project

import numpy as np
import matplotlib.pyplot as plt
import tbplas as tb
from mos2class import mcell, mmetric, mxtohopvec

import ast
#-----------------------------------------------------
#VARIABLES (IMPORTANT):

path_runs = "../../highiterations_newk/"
path_output = "../../importantresults/"

ID = 2111 #IDs of Runs to evaluate
E_min = 0.1 #Must be same as in pca_graddesc.py



#-----------------------------------------------------
k_points = np.array([
        [0.0, 0.0, 0.0],    # Gamma
        [1./2, 0.0, 0.0],   # M
        [1./3, 1./3, 0.0],  # K
        [0.0, 0.0, 0.0],    # Gamma
    ])

k_path, k_idx = tb.gen_kpath(k_points, [40, 40, 40])
k_path_efficient, k_idx_efficient = tb.gen_kpath(k_points, [15, 15, 15]) #kpath with less points for faster calculation
k_label = ["G", "M", "K", "G"]

dset = [0.0, 0.24, 0.32]

widthvector = [0.5,1.0,0.4]

colors = ["red", "black", "grey", "#984ea3", "#ff7f00", "#ffff33", "#a65628"]
#-----------------------------------------------------
#MAIN:

#IDEAL CELL:
ideal_cell = mcell("Ideal",E_min)
idealhops = ideal_cell.mhoppings
#ideal_bandss = ideal_cell.mcalcbands()



#-----------------------------------------------------
#Main:
#-----------------------------------------------------
def safebandstructure(mcellvector,filename):
    import matplotlib.pyplot as plt

    n= len(mcellvector)
    cellvector = [mcell.mprimcell for mcell in mcellvector]

    k_len_vector=[]
    bands_vector=[]
    for i in range(n):
        k_len, bands = tb.calc_bands(cellvector[i],k_path,echo_details=False)
        k_len_vector.append(k_len)
        bands_vector.append(bands)

    plt.figure(dpi = 200)
    #Plotting the Bandstructure
    for i in range(n):
        color = colors[i]
        num_bands = bands_vector[i].shape[1]
        width = widthvector[i]
        for j in range(num_bands):
            if(j==0):
                label = str(mcellvector[i].mname) + " , " + str(mcellvector[i].mnhoppings) +" hoppings"
                plt.plot(k_len_vector[i], bands_vector[i][:, j], color=color, linewidth=width, label = label)
            else:
                plt.plot(k_len_vector[i], bands_vector[i][:, j], color=color, linewidth=width)

        for idx in k_idx:
            plt.axvline(k_len_vector[i][idx], color='k', linewidth=1.0)
        plt.xlim((0, np.amax(k_len_vector[i])))
        plt.xticks(k_len_vector[i][k_idx], k_label)
    plt.xlabel("k (1/nm)")
    plt.ylabel("Energy (eV)")
    #plt.tight_layout()
    plt.legend(
    loc = 'upper center',
    bbox_to_anchor = (0.5,-0.18),
    ncol =1
    )
    plt.subplots_adjust(bottom=0.30)
    pngname = filename + ".png"
    #plt.tight_layout()
    plt.savefig(pngname)
    #plt.show()
    plt.close()


#Reading x from file
name = path_runs + "graddesc_finalx_run" + str(ID) + ".txt"
finalxfile = open(name,"r")
content = finalxfile.read()
content = content.replace("np.float64", "")
x = ast.literal_eval(content)

#Calculating hopvec:
hopvec = mxtohopvec(x,idealhops)
hopvec = [hop for hop in hopvec if abs(hop[3])>E_min]

#Build the final cell:
reducedhopcell = mcell("hoppings reduced by NGD", E_min)
reducedhopcell.mchangehops_tohopvec(hopvec)
N = reducedhopcell.mnhoppings

#Build the "true" cell:
truecell = mcell("original hoppings", 0)
total = truecell.mnhoppings
true_bands = truecell.mcalcbands()


#Plot Bandstructure Comparison with "true" cell
cellvector=[reducedhopcell,truecell]
filename= path_output + "graddesc_bands_run"+str(ID)
title = "band structures of MoS2 cell with all hoppings \n vs MoS2 cell with reduced hoppings, \n using Nesterov-Gradient-Descent"
safebandstructure(cellvector,filename)

#Plot Bandstructure Comparison with reduced and energyreduced cell
reducedenergycell = mcell("hoppings reduced in ascending order of Energy",0)
allhops_array = reducedenergycell.mhoppingtable()[0]
allhops_array = allhops_array[allhops_array[:,6].argsort()[::-1]]
filtered_array = allhops_array[:N]
reducedenergycell.mchangehops_toarr(filtered_array)

cellvector=[reducedhopcell,truecell,reducedenergycell]
filename=path_output + "graddesc_bands_redE_run"+str(ID)
title = "band structures of MoS2 cell with all hoppings \n vs MoS2 cell with Nesterov-GD-reduced hoppings \n vs MoS2 cell with reduced hoppings in order of Energy"
safebandstructure(cellvector,filename)


print("Safed bandstructures succesfully")
