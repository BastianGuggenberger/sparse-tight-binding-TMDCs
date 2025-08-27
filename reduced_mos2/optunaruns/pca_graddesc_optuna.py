#Nesterov Gradient Descent, runs through optuna

#Important for Parallelization:
import os
os.environ["OMP_NUM_THREADS"]       = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"]      = "1"
os.environ["NUMEXPR_NUM_THREADS"]  = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import math
import numpy as np
import time
import optuna
import pandas as pd

#import mos2class
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from Projekte.reduced_mos2.resources.mos2class.mos2class import mcell, safebandstructure, metric, xtohopvec

#CONSTANTS:
run = 27 #Increase by 1 for every run !!!!
path = "optunaresults/" #path where the results will be stored
storage_URL = "sqlite:///optuna_mos2.db"  #URL for the parallelization


def objective(trial):
#-----fixed:

    E_min = 0.2 #Must be same as in results.py

    lambda_0 = 2.3 #Weight of the metric term in the EF
    lambda_1 = 17 #Weight of the sqrt term in the EF
    lambda_2 = 0.05 #Weight of the power 6 term in the EF

    #kappa = 1/1700 #Speed of Gradient Descent
    deltax = 0.001 #deltax in the partial derivative
    iterations = 400 #Iterations of Gradient descent
    #gamma = 0.3 #Factor gamma in the Nesterov acceleration

#-----free:

    kappa = trial.suggest_float("kappa",0.0001,0.002)
    gamma = trial.suggest_float("gamma",0,0.6)


#----------------------
    #IDEAL CELL:
    ideal_cell = mcell("Ideal",E_min)
    idealhops = ideal_cell.mhoppings
    ideal_k_lenn, ideal_bandss = ideal_cell.calcbands()


    #-----------------------------------------------------
    #Functions:
    #-----------------------------------------------------

    def EF(x,currentcell):
        
        N = sum(math.sqrt(abs(xi)) for xi in x)
        H = sum(xi ** 6 for xi in x)
        
        ef = lambda_0 *metric(currentcell,ideal_k_lenn,ideal_bandss) + lambda_1*N + lambda_2*H
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

    x = [1 for hop in idealhops]
    v = [0 for hop in idealhops] #necessary for nesterov gd
    currentcell = mcell("currentcell",E_min)

    for wdh in range(iterations):
        #start = time.time()

        #1.Step Nesterov: x_tilde(t) = x(t) + gamma*v(t-1):
        y = x.copy()
        y = [y[i] + gamma * v[i] for i in range(len(y))]
        
        #2.Step Nesterov: v(t) = gamma*v(t-1) - kappa*d(EF)/dx (x_tilde(t))
        currenthopvec = xtohopvec(y,idealhops)
        currentcell.changehops_tohopvec(currenthopvec)
        Efx = EF(y,currentcell)

        for i in range(len(y)):
            partial = part_deriv_EF(y,i,Efx,currentcell)
            v[i] = gamma * v[i] - kappa * partial
            #3.Step Nesterov: x(t+1) = x(t) + v(t)
            x[i] = x[i] + v[i]

        resultcell = mcell("result",E_min)
        resultcell.changehops_tohopvec(xtohopvec(x,idealhops))

        #Important results:
        N = sum(1 for i in range(len(y)) if abs(idealhops[i][3]*y[i]) <= 0.2)
        try:
            m = metric(resultcell,ideal_k_lenn,ideal_bandss)
        except RuntimeError:
            raise optuna.exceptions.TrialPruned()

        if(N==0):
            loss = m
        else:
            loss = m/N

        #print("Loss = ", loss)
        #Optuna Pruning:
        #if wdh % 10 == 0: #trial report stops parallelization - only report every 10th iteration
            #trial.report(loss, wdh)
            #if trial.should_prune():
                #raise optuna.exceptions.TrialPruned()

    return loss

resultfile = open("results.txt","w+")
resultfile.writelines("Optuna - best params:")

historyfile = open("history.txt","w+")

#For the parallelization:
import sqlite3
conn = sqlite3.connect("optuna_mos2.db")
conn.execute("PRAGMA journal_mode = WAL;")
conn.close()

study = optuna.create_study(
    study_name="pca_graddesc_3",
    storage=storage_URL,
    load_if_exists=True,       # lädt eine bestehende Study, statt sie neu anzulegen
    direction="minimize",
    pruner=optuna.pruners.MedianPruner(
        n_startup_trials=5,
        n_warmup_steps=2
    ),)


if(__name__ == "__main__"):
    study.optimize(objective,n_trials = 200, n_jobs=4)
    df = study.trials_dataframe()
    df.to_csv("optuna_trials.csv", index = False)
    print(study.best_params)
    print("Loss=", study.best_value)

resultfile.close()
historyfile.close()