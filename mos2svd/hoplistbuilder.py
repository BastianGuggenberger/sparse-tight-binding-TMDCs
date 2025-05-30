import tbplas as tb
import sys
import io
import ast

#mhoppinglist(): returns a list of hopping informations for all hoppings in the testcell
def mhoppinglist(testcell):
    hoppinglist = []
    orbsetdict = {}
    count = 0 #Counts hoppings which are also in the lower triangle

    #output = io.StringIO()
    #Versuch:
    old_stdout = sys.stdout
    try:
        output = io.StringIO()
        sys.stdout = output
        testcell.print()
        text = output.getvalue()
    finally:
        sys.stdout = old_stdout

    #sys.stdout = output
    #testcell.print()
    #sys.stdout = sys.__stdout__
    
    #text = output.getvalue()
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

    print("Hoppings not added: ", count)

    return hoppinglist

idealcell = tb.make_mos2_soc()
idealhoppingfile = open("idealhops.txt","w+")
hoplist = mhoppinglist(idealcell)
for i in range(len(hoplist)):
    idealhoppingfile.write(str(hoplist[i]))
    if (i != len(hoplist)-1) :
        idealhoppingfile.write("\n")
idealhoppingfile.close()