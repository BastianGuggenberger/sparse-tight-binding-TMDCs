#Zwei verschiedene Möglichkeiten eine 1D PCA durchzufüren:
#(Ausgehend vom der 4D Iris Dataset)
#Möglichkeit 1: 1D PCA durchführen
#Möglichkeit 2: Zuerst 3D PCA durchführen und anschliesend dieses 3D Dataset zu 1D PCA analysieren
#Theorie: Kein Unterschied

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from sklearn.decomposition import PCA



#Datenbeschaffung

Rohdaten = load_iris().data
Rohdaten_Detailliert = load_iris(as_frame=True)
print(Rohdaten.shape)




#PCAs Durchführen

Daten_1D_direkt = PCA(n_components=1).fit_transform(Rohdaten)
Daten_3D = PCA(n_components=3).fit_transform(Rohdaten)
Daten_1D_indirekt = PCA(n_components=1).fit_transform(Daten_3D)

print(Daten_1D_direkt)



#Vergleich

vergleichsarray = np.empty((0,3))
for i in range(Daten_1D_direkt.shape[0]):
    vergleichsarray = np.vstack((vergleichsarray, [Daten_1D_direkt[i][0], Daten_1D_indirekt[i][0], Daten_1D_direkt[i][0]-Daten_1D_indirekt[i][0]]))

print(vergleichsarray)



#Graphik

fig1 = plt.figure()
plt1 = fig1.add_subplot(121)
plt2 = fig1.add_subplot(122)
DatenA = np.empty((0,2))
DatenB = np.empty((0,2))
for i in range(Daten_1D_direkt.shape[0]):
    DatenA = np.vstack((DatenA,[20,Daten_1D_direkt[i][0]]))
    DatenB = np.vstack((DatenB,[20,Daten_1D_indirekt[i][0]]))

plt1.scatter(
    DatenA[:,0],
    DatenA[:,1],
   c = Rohdaten_Detailliert.target,
    s = 2
)
plt2.scatter(
    DatenB[:,0],
    DatenB[:,1],
    c = Rohdaten_Detailliert.target,
    s = 2
)

plt.show()