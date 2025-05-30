from svdmos2class import msvdcell, metric, safebandstructure
E_min = 0.1

for N in range(2,5):
    idealcell = msvdcell("Ideal Cell",E_min)
    svdcell = msvdcell("SVD based hoppings reduced Cell",E_min)

    print("Original Hoppings: ", svdcell.mnhoppings)
    svdcell.mtruncate(N)
    print("New Hoppings: ",svdcell.mnhoppings)

    cellvec = [idealcell,svdcell]
    safebandstructure(cellvec,"comparison_N=" + str(N), "Comparison of SVD reduced vs Ideal Hoppings, N = " + str(N))