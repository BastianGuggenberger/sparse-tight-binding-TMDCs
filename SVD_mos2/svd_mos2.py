#Compares the bandstructures of SVD based reduced MoS2 cell
#iterates for different n_components

from mos2class  import mcell, msafebandstructure

#-----------------------------------------------------------------------
#VARIABLES: (Important!)
oldsvdoutputpath = "results/oldsvd/"
newsvdoutputpath = "results/newsvd/"
#-----------------------------------------------------------------------
E_min = 0.1
idealcell = mcell("Ideal Cell",E_min)

#Old SVD:
for N in range(1,5):
    svdcell = mcell("old SVD",E_min)
    svdcell.msvdtruncate(N)
    path = oldsvdoutputpath + "graddesc_finalhopvec_old"+str(N)+".txt"
    finalhopfile = open(path, "w+")
    finalhopfile.write(str(svdcell.mhoppings.copy()))
    finalhopfile.close()


#New SVD:
for N in range(1,23):
    svdcell = mcell("new SVD",E_min)
    svdcell.mseperatedsvdtruncate(N)
    path = newsvdoutputpath + "graddesc_finalhopvec_new"+str(N)+".txt"
    finalhopfile = open(path, "w+")
    finalhopfile.write(str(svdcell.mhoppings.copy()))
    finalhopfile.close()