#Analyses and plots neighbours degrees of the original Roldán et al tight binding model

import matplotlib.pyplot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mos2class import mcell, msafebandstructure, mmetric, mxtohopvec

#-----------------------------------------------------
path_output = "pngs/degrees.png"
E_min = 0.1

#-----------------------------------------------------
#Main:
#-----------------------------------------------------

cell = mcell("test", E_min)
orders = cell.mget_neighbours_orders()
print(orders)


labels = ["1. Order", "2. Order"]
colors = ["lightcoral","cornflowerblue"]

fig, ax = plt.subplots(figsize=(6,4))
ax.bar(labels, [orders[1],orders[2]], width = 0.6 , color = colors)
ax.set_ylabel("Number of hoppings")
#ax.set_xlim(-0.2,1.2)
plt.savefig(path_output)