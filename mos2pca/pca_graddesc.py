#Script for Nesterov Gradient Descent based PCA of MoS2 Hopping Terms
#Based on the Class "mos2class.py"
#Based on tbplas method make_mos2_soc, which is an implementation of:
    # R Roldán et al 2014 2D Mater. 1 034003
    # https://www.tbplas.net/_api/tbplas.make_mos2_soc.html

import math
import numpy as np
from mos2class import mcell, mmetric, mxtohopvec
import time

#CONSTANTS:

lambda_1vec = []
runvec = [] #Increase by 1 for every run !!!! #Declared
resultpath = "results/highiterations_2.0/" #path where the results will be stored

E_min = 0.1 #Must be same as in results.py

#Weights: Change for different m and N results
lambda_0 = 4.0 #Weight of the metric term in the EF
#lambda_1 = 50.0 #Weight of the sqrt term in the EF

#Hyperparameters (keep constant):
startingx = 1
lambda_2 = 0.05 #Weight of the power 6 term in the EF
kappa = 1/1700 #Speed of Gradient Descent
deltax = 0.05 #deltax in the partial derivative
iterations = 1200 #Iterations of Gradient descent
gamma = 0.3 #Factor gamma in the Nesterov acceleration


#IDEAL CELL:
ideal_cell = mcell("Ideal",E_min)
idealhops = ideal_cell.mhoppings.copy()
ideal_bands = ideal_cell.mcalcbands()
ideal_bands_efficient = ideal_cell.mcalcbands(efficient=True)


#-----------------------------------------------------
#Functions:
#-----------------------------------------------------

#errorfunction that should be minimized by graddesc
# EF = lambda_0 * metric + lambda_1 * sqrt(|x|) + lambda_2 * x^6
def EF(x,currentcell):
    
    N = sum(math.sqrt(abs(xi)) for xi in x)
    H = sum(xi ** 6 for xi in x)
    
    metric = mmetric(currentcell,ideal_bands_efficient,efficient=True)
    

    ef = lambda_0 * metric + lambda_1*N + lambda_2*H

    return ef

#partial derivative of the errorfuncion, in the direction of hop i
#calculated numerically
def part_deriv_EF(x,i,Efx,currentcell):
    xnew = x.copy()
    xnew[i] += deltax
    old_E = currentcell.mhoppings[i][3]
    currentcell.mchangehop_energy(i,xnew[i]*idealhops[i][3])
    errorfunction = EF(xnew,currentcell)
    diff = (errorfunction - Efx) / deltax
    currentcell.mchangehop_energy(i,old_E)

    return diff



#-----------------------------------------------------
#Main:
#-----------------------------------------------------

def nesterovgd(run,lambda_1,printinfos=True):
    #File management
    end = str(run) + ".txt"
    historyfile = open(resultpath + "graddesc_history_run"+end, 'w+')
    xfile = open(resultpath + "graddesc_x_run"+end, 'w+')
    finalxfile = open(resultpath + "graddesc_finalx_run"+end, 'w+')
    paramsfile = open(resultpath + "graddesc_params_run"+end, 'w+')

    #Note parameters
    paramsfile.writelines("run = " + str(run) + "\n")
    paramsfile.writelines("E_min = " + str(E_min) + "\n")
    paramsfile.writelines("lambda_0 = " + str(lambda_0) + "\n")
    paramsfile.writelines("lambda_1 = " + str(lambda_1) + "\n")
    paramsfile.writelines("lambda_2 = " + str(lambda_2) + "\n")
    paramsfile.writelines("kappa = " + str(kappa) + "\n")
    paramsfile.writelines("deltax = " + str(deltax) + "\n")
    paramsfile.writelines("iterations = " + str(iterations) + "\n")
    paramsfile.writelines("gamma = " + str(gamma) + "\n")
    paramsfile.close()


    #NESTEROV GRADIENT DESCENT:
    x = [startingx for hop in idealhops] #weight vector x, gives the relative weight for the hopping energy of each hopping
    v = [0 for hop in idealhops] #velocity vector v in nesterov gradient descent
    currentcell = mcell("currentcell",E_min)

    orderstotal = currentcell.mget_neighbours_orders()
    orderstotal_sum = sum(orderstotal)

    for wdh in range(iterations):
        #start = time.time()

        #1.Step Nesterov: x_tilde(t) = x(t) + gamma*v(t-1):
        y = x.copy()
        y = [y[i] + gamma * v[i] for i in range(len(y))]
        
        #2.Step Nesterov: v(t) = gamma*v(t-1) - kappa*d(EF)/dx (x_tilde(t))
        currenthopvec = mxtohopvec(y,idealhops.copy())
        currentcell.mchangehops_tohopvec(currenthopvec)
        Efx = EF(y,currentcell)

        for i in range(len(y)):
            partial = part_deriv_EF(y,i,Efx,currentcell)

            v[i] = gamma * v[i] - kappa * partial
            #3.Step Nesterov: x(t+1) = x(t) + v(t)
            x[i] = x[i] + v[i]


        resultcell = mcell("result",E_min)
        resultcell.mchangehops_tohopvec(mxtohopvec(x,idealhops.copy()))
        N_1 = sum(1 for i in range(len(y)) if ((abs(idealhops[i][3]*y[i]) <= E_min) and (currentcell.mget_neighbours_order(idealhops.copy()[i])==1)))
        N_2 = sum(1 for i in range(len(y)) if ((abs(idealhops[i][3]*y[i]) <= E_min) and (currentcell.mget_neighbours_order(idealhops.copy()[i])==2)))
        #end = time.time()

        #Print iteration results:
        string = ""
        string += "n = " + str(wdh) + "\n"
        string += "N(x_i(1) ~ 0) = " + str(N_1) + " / " + str(orderstotal[1]) +"\n"
        string += "N(x_i(2) ~ 0) = " + str(N_2) + " / " + str(orderstotal[2]) +"\n"
        string += "metric=" + str(mmetric(resultcell,ideal_bands)) + "\n" + "\n" #Here the less efficient, more accurate metric is calculated
        
        if(printinfos==True):
            print(string)
        historyfile.writelines(string)
        xfile.writelines(str(x))


    finalxfile.writelines(str(x))
    historyfile.close()
    xfile.close()
    finalxfile.close()

runlambdapairs = zip(runvec,lambda_1vec)
for pair in runlambdapairs:
    run = pair[0]
    lambda_1 = pair[1]
    print("Run " + str(run) + " in progress.")
    nesterovgd(run,lambda_1,False)
    print("Run " + str(run) + " done.")
