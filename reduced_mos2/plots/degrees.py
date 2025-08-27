#Analying Neighbours degrees of the blueprint tight binding model

import matplotlib.pyplot
import ast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#import mos2class
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from Projekte.reduced_mos2.resources.mos2class.mos2class import mcell, msafebandstructure, mmetric, mxtohopvec

#-----------------------------------------------------
#path_output = "importantresults_2.0/"
E_min = 0.1

#-----------------------------------------------------
#Main:
#-----------------------------------------------------

cell = mcell("test", E_min)
orders = cell.mget_neighbours_orders()
print(orders)


labels = ["1. Order", "2. Order"]
colors = ["lightcoral","cornflowerblue"]

fig, ax = plt.subplots(figsize=(2,4))
ax.bar(labels, [orders[1],orders[2]], width = 0.6 , color = colors)
ax.set_ylabel("Number of hoppings")
#ax.set_xlim(-0.2,1.2)
plt.show()