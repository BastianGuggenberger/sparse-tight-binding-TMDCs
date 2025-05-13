import math
import numpy as np
from mos2class import mcell, metric
import time

#CONSTANTS:
E_min = 0.2

lambda_0 = 1.6 #Weight of the metric term in the EF
lambda_1 = 7 #Weight of the sqrt term in the EF
lambda_2 = 4 #Weight of the power 6 term in the EF

kappa = 1/400 #Speed of Gradient Descent
deltax = 0.1 #deltax in the partial derivative
iterations = 150 #Iterations of Gradient descent

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

def EF(x,currentcell):
    
    N = sum(math.sqrt(abs(xi)) for xi in x)
    H = sum(xi ** 6 for xi in x)
    
    ef = lambda_0 *metric(ideal_cell,currentcell) + lambda_1*N + lambda_2*H
    return ef

def part_deriv_EF(x,i,Efx,currentcell):
    xnew = x.copy()
    xnew[i] += deltax
    old_E = currentcell.mhoppings[i][3]
    currentcell.mchangehop_energy(i,xnew[i]*idealhops[i][3])
    diff = (EF(xnew,currentcell) - Efx) / deltax
    currentcell.mchangehop_energy(i,old_E)

    return diff



#-----------------------------------------------------
#Main:
#-----------------------------------------------------
historyfile = open("graddesc_history.txt", 'w+')
xfile = open("graddesc_x.txt", 'w+')


x = [1 for hop in idealhops]
currentcell = mcell("currentcell",E_min)

for wdh in range(iterations):
    y = x.copy()
    currenthopvec = xtohopvec(y)
    currentcell.changehops_tohopvec(currenthopvec)
    Efx = EF(y,currentcell)

    for i in range(len(y)):
        partial = part_deriv_EF(y,i,Efx,currentcell)
        x[i] -= kappa * partial

    resultcell = mcell("result",E_min)
    resultcell.changehops_tohopvec(xtohopvec(x))
    N = sum(1 for i in range(len(y)) if abs(idealhops[i][3]*y[i]) <= 0.2)
    """
    print("n = ", wdh)
    print("N(x_i ~ 0) = ",N)
    print("metric=" ,metric(resultcell,ideal_cell), "\n")
    """
    string = ""
    string += "n = " + str(wdh) + "\n"
    string += "N(x_i ~ 0) = " + str(N) + "\n"
    string += "metric=" + str(metric(resultcell,ideal_cell)) + "\n" + "\n"
    print(string)
    historyfile.writelines(string)
    xfile.writelines(str(x))

historyfile.close()
xfile.close()
