import matplotlib.pyplot as plt

#import mos2class
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from Projekte.reduced_mos2.resources.mos2class.mos2class import mcell

import ast
#-----------------------------------------------------
#Variables:

historypath = "highiterations_2.0/graddesc_history_run"
outputpath = "lowiterations_2.0/mvsNanalysis/"

typefile = "seplevels"
runvec = range(1000,1034)
rundifference = 100 #how much higher should the output ids be

iterations = 399

E_min = 0.1
#-----------------------------------------------------
#MAIN:

truecell = mcell("True Cell", 0)
total = truecell.mnhoppings

idealcell = mcell("Ideal", E_min)
idealhops = idealcell.mnhoppings

for ID in runvec:

    historyfile = open(historypath + str(ID) + ".txt","r")
    outputfile = open(outputpath + "mvsN_run" + str(ID+rundifference)+ ".txt","w+")
    #-----------------------------------------------------
    rightline = False
    N_red = None
    m = None
    if(typefile == "alllevels"):
        lineid = 1
        for line in historyfile:
            if(rightline == True):
                if (lineid == 1):
                    x = line.split("=")
                    N_red = int(x[1])
                    lineid = 2
                elif (lineid ==2):
                    x = line.split("=")
                    m = float(x[1])
                    break
            else:
                if(line[0]=="n"):
                    x = line.split("=")
                    n = int(x[1])
                    if(n==iterations):
                        rightline=True
    elif (typefile == "seplevels"):
        lineid = 1
        for line in historyfile:
            if(rightline == True):
                if (lineid == 1):
                    x = line.split("=")
                    y = x[1].split("/")
                    N_red = int(y[0])
                    lineid = 2
                elif (lineid == 2):
                    x = line.split("=")
                    y = x[1].split("/")
                    N_red += int(y[0])
                    lineid = 3
                elif (lineid ==3):
                    x = line.split("=")
                    m = float(x[1])
                    break
            else:
                if(line[0]=="n"):
                    x = line.split("=")
                    n = int(x[1])
                    if(n==iterations):
                        rightline=True

    outputfile.write("["+str(float((idealhops-N_red))/total)+"," +str(m)+"]")
    historyfile.close()
    outputfile.close()
    print("Information of run "+ str(ID)+ " safed.")

