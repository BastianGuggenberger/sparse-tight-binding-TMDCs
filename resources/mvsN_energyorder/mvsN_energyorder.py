#Building Primitive Cells with Hoppings in certain Energy intervals
#Working with the mcell class
#Based on Custom Mos2 Class

from mos2class.mos2class import mcell, msafebandstructure, mmetric, mxtohopvec
import tbplas as tb
import numpy as np
import matplotlib.pyplot as plt
#-----------------------------------------------------

#------------------------------------------------------
#1.: Plotting e_min vs remaining hoppings
truecell = mcell("truecell",0)
truebands = truecell.mcalcbands()

minhopa = 0.0
minhopb = 1.6201
total = truecell.mnhoppings
nvec = []
metrics = []
minhops = np.linspace(minhopa,minhopb,100)

for minhop in minhops:
    testcell = mcell("testcell",minhop)
    n = testcell.mnhoppings
    metric = mmetric(testcell,truebands)
    nvec.append(n/total)
    metrics.append(metric)
    output = "totalhops: " + str(testcell.mnhoppings) + ",  relhops: " + str(testcell.mnhoppings/total) + ",  E_min: " + str(minhop) + ",  metric: " + str(mmetric(testcell,truebands))
    print(output)

#Save the data of the figure
mvsNfile = open("mvsN_energyorder.txt",'w+')
mvsNfile.writelines(str(nvec)+"\n")
mvsNfile.writelines(str(metrics))
mvsNfile.close()

#Plotting:
plt.gca().invert_xaxis()
plt.axhline(y=0,color='grey')
plt.plot(nvec, metrics, color = "Red", label="Hoppings reduced in the order of increasing energies")
plt.ylim(top=500,bottom=0) #Important!
plt.xlabel("Number N of Hoppings / Total Number of Hoppings")
plt.ylabel("Error Metric m in eV")
plt.title("Error Metric m for varying number N of hoppings. \n lambda_0 = ")
plt.savefig("mvsnenergyorder.png")
plt.clf()