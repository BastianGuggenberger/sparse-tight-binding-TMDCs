#Script for analyzing and plotting Results of Test runs of the Nesterov Gradient Descent PCA
#Based on the script "pca_graddesc.py"
import numpy as np
import matplotlib.pyplot as plt
import tbplas as tb

#import mos2class
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from Projekte.reduced_mos2.resources.mos2class.mos2class import mcell, mmetric, mxtohopvec

import ast
#-----------------------------------------------------
#VARIABLES (IMPORTANT):

path_output = "pngs/bandgapillustration"

ID = 1007 #IDs of Runs to evaluate
E_min = 0.1 #Must be same as in pca_graddesc.py



#-----------------------------------------------------
k_points = np.array([
        [0.0, 0.0, 0.0],    # Gamma
        [1./2, 0.0, 0.0],   # M
        [2./3, 1./3, 0.0],  # K
        [0.0, 0.0, 0.0],    # Gamma
    ])

k_path, k_idx = tb.gen_kpath(k_points, [40, 40, 40])
k_path_efficient, k_idx_efficient = tb.gen_kpath(k_points, [15, 15, 15]) #kpath with less points for faster calculation
k_label = ["G", "M", "K", "G"]

colors = [
    "#e41a1c", 
    "#377eb8", 
    "#4daf4a", 
    "#984ea3", 
    "#ff7f00",  
    "#ffff33", 
    "#a65628",  
    "#f781bf", #Diese
    "#999999", 
    "#17becf", #Diese
    "#e41a1c", 
    "#377eb8", 
    "#4daf4a", 
    "#984ea3", 
    "#ff7f00",  
    "#ffff33", 
    "#a65628",  
    "#f781bf", 
    "#999999", 
    "#17becf",
    "#e41a1c", 
    "#377eb8", 
    "#4daf4a", 
    "#984ea3", 
    "#ff7f00",  
    "#ffff33", 
    "#a65628",  
    "#f781bf", 
    "#999999", 
    "#17becf",
]

#-----------------------------------------------------
#Main:
#-----------------------------------------------------
def safebandstructure(cell,filename):
    import matplotlib.pyplot as plt
    

    k_len, bands = cell.mprimcell.calc_bands(k_path,echo_details=False)


    plt.figure(dpi = 200)
    #Plotting the Bandstructure
    num_bands = bands.shape[1]
    for j in range(num_bands):
        if(j>9):
            plt.plot(k_len, bands[:, j], color="blue", linewidth="1.0")
        elif(j<7):
            plt.plot(k_len, bands[:, j], color="green", linewidth="1.0")
        elif(j==9):
            plt.plot(k_len, bands[:, j], color="blue", linewidth="1.0", label = "MoS2 valence bands")
        elif(j==7):
            plt.plot(k_len, bands[:, j], color="green", linewidth="1.0", label = "MoS2 conduction bands")

    #7,9
    for idx in k_idx:
        plt.axvline(k_len[idx], color='k', linewidth=1.0)
    plt.axhline(-2.0, color = "k", linewidth = 1.0, label = "E = -2.0 eV / E = 2.0 eV")
    plt.axhline(2.0, color = "k", linewidth = 1.0)

    plt.xlim((0, np.amax(k_len)))
    plt.xticks(k_len[k_idx], k_label)
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


#Build the "true" cell:
truecell = mcell("original hoppings", 0)
total = truecell.mnhoppings
true_bands = truecell.mcalcbands()
safebandstructure(truecell,path_output)

print("Plotted BS illustration succesfully")
