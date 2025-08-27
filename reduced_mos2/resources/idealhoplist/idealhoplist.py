#idealhoplist.py creates the file "idealhoplist.txt" containing hoppings as rows in python-list format.  
#The stored hoppings are all hoppings from the  R Roldán et al paper without their imaginary part.  
#The TBPLAS-implementation of the R Roldán model "make_mos2_soc" is used for that, however this implementation contains some hoppings twice, which is also fixed here.  
#The resulting hoplist ist used each time a new "mcell" is constructed.  

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
idealhoppingfile = open("idealhoplist.txt","w+")
hoplist = mhoppinglist(idealcell)
for i in range(len(hoplist)):
    idealhoppingfile.write(str(hoplist[i]))
    if (i != len(hoplist)-1) :
        idealhoppingfile.write("\n")
idealhoppingfile.close()

print("Succesfully stored hoppings in 'idealhoplist.txt'")
