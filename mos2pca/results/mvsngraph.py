#mvsngraph
import numpy as np
import matplotlib.pyplot as plt
import ast

path = "Nesterovruns/mvsNanalysis/"

Nvec = []
mvec = []
Nvsm=[]
for ID in range(12,26):
    if (ID == 14):
        continue
    name = path + "mvsN_run" + str(ID) + ".txt"
    finalxfile = open(name,"r")
    content = finalxfile.read()
    x = ast.literal_eval(content)
    finalxfile.close()
    Nvsm.append(x)
    Nvec.append(x[0])
    mvec.append(x[1])

#Sorting
Nvsm = np.array(Nvsm)
Nvsm = Nvsm[Nvsm[:,0].argsort()]

#Get N vs M Information of Energyreduced Cell:
name = path + "mvsN_fromlowEtohighE.txt"
finalxfile = open(name,"r")
content = []
for line in finalxfile:
    content.append(line)
finalxfile.close()
ns = ast.literal_eval(content[0])
metrics = ast.literal_eval(content[1])

#Plotting:
plt.gca().invert_xaxis()
plt.axhline(y=0,color='grey')
plt.plot(ns, metrics, color = "red", label="Hoppings reduced in the order of increasing energies")
plt.plot(Nvsm[:,0], Nvsm[:,1], color = "green", label="Hoppings reduced with Nesterov Gradient Descent")
plt.ylim(top=3000) #Important!
plt.xlabel("Number N of Hoppings / Total Number of Hoppings")
plt.ylabel("Error Metric m")
plt.legend()
plt.title("Error Metric m for different number N of hoppings.")
plt.savefig(path + "Nesterov_Nvsm"+".png")
plt.clf()