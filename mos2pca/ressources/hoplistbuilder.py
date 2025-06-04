#Based on tbplas method make_mos2_soc, which is an implementation of:
    # R Roldán et al 2014 2D Mater. 1 034003
    # https://www.tbplas.net/_api/tbplas.make_mos2_soc.html

#hoplistbuilder.py stores all hoppings of the R Roldán et al paper without their imaginary part.
#some hoppings of the make_mos2_soc implementation are added twice, which is also fixed here.

import tbplas as tb
import sys
import io
import ast

#mhoppinglist(): returns a list of all hoppings in the implementation without imaginary part.
#formatted in a way that is useful for the mcell class
def mhoppinglist(cell):
    hoppinglist = []
    count = 0 #Counts hoppings which are also in the lower triangle of the hamiltonian matrix

    old_stdout = sys.stdout
    try:
        output = io.StringIO()
        sys.stdout = output
        cell.print()
        text = output.getvalue()
    finally:
        sys.stdout = old_stdout

    text = text.split("terms:",1)[1]
    for line in text.splitlines():
        if(line == ""):
            continue
        for char in ("(",")",","):
            line = line.replace(char, "%")
        linevec = line.split("%")
        linevec = [element for element in linevec if element.strip()]
        rn = [int(linevec[0]),int(linevec[1]),int(linevec[2])]
        orb_i = int(linevec[3])
        orb_j = int(linevec[4])

        #Energy: (Deletes imaginary parts)
        energystring = linevec[5]
        energycomplex = ast.literal_eval(energystring.strip())
        energyreal = energycomplex.real

        #Check for unnecesaary hoppings (conjugate already in hamiltonian)
        hopunique = True
        for hop in hoppinglist:
            if(hop[0][0]==rn[0] and hop[0][1]==rn[1] and hop[0][2]==rn[2]):
                if((hop[1]==orb_i and hop[2]==orb_j) or (hop[1]==orb_j and hop[2]==orb_i)):
                    hopunique = False

        if(hopunique == True):
            hoppinglist.append([rn,orb_i,orb_j,energyreal])
        else:
            count+=1

    print("Non-unique hoppings (not added): ", count)

    return hoppinglist


#Write down hoppings in hoppingslist
idealcell = tb.make_mos2_soc()
idealhoppingfile = open("idealhops.txt","w+")
hoplist = mhoppinglist(idealcell)
for i in range(len(hoplist)):
    idealhoppingfile.write(str(hoplist[i]))
    if (i != len(hoplist)-1) :
        idealhoppingfile.write("\n")
idealhoppingfile.close()

print("Succesfully stored hoppings in 'idealhops.txt'")
