#Analying Neighbours degrees of multiple hopping reduced cells

import matplotlib.pyplot
import ast
import numpy as np
import pandas as pd

#import mos2class
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from mos2class import mcell, msafebandstructure, mmetric, mxtohopvec

#-----------------------------------------------------

path_graddesc = "finalruns/graddesc_finalx_run40.txt"
path_output = "finalruns/formated_results/degree_analysis/"

E_min = 0.1 #Must be same as in the graddesc and svd runs

#-----------------------------------------------------
#Main:
#-----------------------------------------------------

truecell = mcell("all hoppings", 0)
idealcell = mcell("ideal", E_min)
idealhops = idealcell.mhoppings.copy()

#Build graddesccell

#Reading x from file
finalxfile = open(path_graddesc,"r")
content = finalxfile.read()
x = ast.literal_eval(content)

#Calculating hopvec:
hopvec = mxtohopvec(x,idealhops)
hopvec = [hop for hop in hopvec if abs(hop[3])>E_min]

#Build the final graddesc cell:
graddeschopcell = mcell("hoppings reduced by Nesterov grad-desc", E_min)
graddeschopcell.mchangehops_tohopvec(hopvec)
N_graddesc = graddeschopcell.mnhoppings



#Get the table:
table = [truecell.mget_neighbours_orders(),idealcell.mget_neighbours_orders(),graddeschopcell.mget_neighbours_orders()]

deg0 = [table[0][0],table[1][0],table[2][0]]
deg1 = [table[0][1],table[1][1],table[2][1]]
deg2 = [table[0][2],table[1][2],table[2][2]]

data={"Hoppings Type": ["All Hoppings","E_min reduced","GradDesc reduced - "+str(N_graddesc)+ " hoppings"],"Degree 0": deg0, "Degree 1": deg1,"Degree 2": deg2,"Degree 3":[0 for x in table]}
df = pd.DataFrame(data)
output = path_output + "Neighbours_Degrees.csv"
df.to_csv(output)

finalxfile.close
