import math
import numpy as np
from mos2class import mcell, metric
import time

#CONSTANTS:
run = 7 #Increase by 1 for every run

E_min = 0.2 #Must be same as in results.py

lambda_0 = 1.0 #Weight of the metric term in the EF
lambda_1 = 25 #Weight of the sqrt term in the EF
lambda_2 = 0.2 #Weight of the power 6 term in the EF

kappa = 1/1000 #Speed of Gradient Descent
deltax = 0.1 #deltax in the partial derivative
iterations = 200 #Iterations of Gradient descent

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
end = str(run) + ".txt"
historyfile = open("graddesc_history_run"+end, 'w+')
xfile = open("graddesc_x_run"+end, 'w+')
finalxfile = open("graddesc_finalx_run"+end, 'w+')
paramsfile = open("graddesc_params_run"+end, 'w+')

#Note parameters
paramsfile.writelines("run = " + str(run) + "\n")
paramsfile.writelines("E_min = " + str(E_min) + "\n")
paramsfile.writelines("lambda_0 = " + str(lambda_0) + "\n")
paramsfile.writelines("lambda_1 = " + str(lambda_1) + "\n")
paramsfile.writelines("lambda_2 = " + str(lambda_2) + "\n")
paramsfile.writelines("kappa = " + str(kappa) + "\n")
paramsfile.writelines("deltax = " + str(deltax) + "\n")
paramsfile.writelines("iterations = " + str(iterations) + "\n")
paramsfile.close()

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

    string = ""
    string += "n = " + str(wdh) + "\n"
    string += "N(x_i ~ 0) = " + str(N) + "\n"
    string += "metric=" + str(metric(resultcell,ideal_cell)) + "\n" + "\n"
    print(string)
    historyfile.writelines(string)
    xfile.writelines(str(x))

finalxfile.writelines(str(x))
historyfile.close()
xfile.close()
finalxfile.close()
