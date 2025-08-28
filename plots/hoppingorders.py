#plots how 1st and 2nd order hoppings look like in the primitive cell

import tbplas as tb
from mos2class import mcell, mxtohopvec

#-----------------------------------------------------
E_min = 0.1
#-----------------------------------------------------

idealcell = mcell("ideal",E_min)
idealhops = idealcell.mhoppings.copy()
order_1hops = []
order_2hops = []
for hop in idealhops:
    if(idealcell.mget_neighbours_order(hop)==1):
        order_1hops.append(hop)
    elif(idealcell.mget_neighbours_order(hop)==2):
        order_2hops.append(hop)

idealcell.mchangehops_tohopvec(order_1hops)
idealcell.mprimcell.plot(fig_name = "pngs/order1hops.png",hop_color = "red")


idealcell.mchangehops_tohopvec(order_2hops)
idealcell.mprimcell.plot(fig_name = "pngs/order2hops.png",hop_color = "blue")