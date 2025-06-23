#Script for analyzing the final x vector and espacially its orders
#The code is based on parts of the code in "results.py" (plotting of the final x vector)

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

path_runs = "highiterations/"
path_output = "importantresults/"
ID = 61
E_min = 0.1


#-----------------------------------------------------
#MAIN:

#analyze oders of ideal cell:
idealcell = mcell("ideal", E_min)
idealhops = idealcell.mhoppings.copy()
orders = []
for hop in idealhops:
    orders.append(idealcell.mget_neighbours_order(hop))

#Reading x from file
name = path_runs + "graddesc_finalx_run" + str(ID) + ".txt"
finalxfile = open(name,"r")
content = finalxfile.read()
x = ast.literal_eval(content)


#ordering x-vector based on orders:
pairs = list(zip(x,orders))
order_1 = [[xi, flag] for xi, flag in pairs if (flag==1)]
order_2 = [[xi, flag] for xi, flag in pairs if (flag==2)]
pairs_ordered = order_1 + order_2
x_ordered = [xi for xi, flag in pairs_ordered]

#Plot ordered x vector
indices = [i for i in range(len(x))]
plt.figure(dpi = 200)
plt.plot(indices,x_ordered,color = "Grey",linewidth = 0.6)
for i, pair in enumerate(pairs_ordered):
    if (pair[1] == 1): color = "Red"
    elif (pair[1] == 2): color = "Blue"
    else: print("order not recognized error.")
    plt.scatter(i,pair[0],color=color)
plt.xlabel("index i")
plt.ylabel("weight factor x[i]")
#plt.title("Weight vector x")
plt.savefig(path_output + "graddesc_x_ordered"+str(ID)+".png")
plt.clf()