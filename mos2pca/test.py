from mos2class import mcell
E_min = 0.1

zellea = mcell("A",E_min)
zelleb = mcell("B",E_min)

x = zellea.mhoppings.copy()
x[3][3] = 45
zellea.changehops_tohopvec(x)
print(zellea.mhoppings[3])

zelleb.mchangehop_energy(3,45)
print(zelleb.mhoppings[3])