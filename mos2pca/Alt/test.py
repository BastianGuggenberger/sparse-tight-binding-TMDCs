#test for mistake with imaginary parts
from mos2class import mcell, safebandstructure

mos2neu = mcell("test",0)
cellvector=[mos2neu]
safebandstructure(cellvector,"test.png","test")