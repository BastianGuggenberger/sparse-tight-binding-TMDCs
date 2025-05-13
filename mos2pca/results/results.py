
import matplotlib.pyplot as plt
import tbplas as tb

#import mos2class
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from mos2class import mcell, safebandstructure

import ast
#-----------------------------------------------------


ID = 3
E_min = 0.2


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
reducedhopcell = mcell("reducedhopcell", E_min)
reducedhopcell.changehops_tohopvec(hopvec)

cellvector=[ideal_cell,reducedhopcell]
filename="graddesc_bands_run"+str(ID)
title = "Bandstructures of ideal mos2cell vs mos2cell with reduced hoppings"
safebandstructure(cellvector,filename,title)
