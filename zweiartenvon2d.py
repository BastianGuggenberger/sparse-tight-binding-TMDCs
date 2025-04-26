#Zwei Arten einer 2-D- PCA:
#1. Variante: 3D PCA von der nur die ersten zwei Komponenten verwendet werden
#2. Variante: 2D PCA
#Der Theorie nach sollte exakt das selbe herauskommen

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from sklearn.decomposition import PCA




#1. DATENBESCHAFFUNG

Rohdaten = load_iris().data
#print(Rohdaten) #Rohdaten ist ein 150x4 numpy array

PCA_3D = PCA(n_components=3).fit_transform(Rohdaten)
PCA_2D = PCA(n_components=2).fit_transform(Rohdaten)





#2. VERGLEICH

vergleich = np.empty((0,4))
for i in range(PCA_3D.shape[0]):
    vergleich= np.vstack((vergleich, [PCA_3D[i][0],PCA_2D[i][0],PCA_3D[i][1],PCA_2D[i][1]]))
print(vergleich)

differenz = np.empty((0,2))
for i in range(PCA_2D.shape[0]):
    differenz = np.vstack((differenz, [vergleich[i][0]-vergleich[i][1], vergleich[i][2]-vergleich[i][3]]))
print(differenz)




#3. GRAPHISCHER VERGLEICH
fig = plt.figure(1, figsize=(8, 6))
plt1 = fig.add_subplot(121, projection="3d")
alledaten = load_iris(as_frame=True)
plt1.scatter(
    PCA_3D[:,0],
    PCA_3D[:,1],
    PCA_3D[:,2],
    s=40,
    c=alledaten.target
)

plt2 = fig.add_subplot(122)
plt2.scatter(
    PCA_2D[:,0],
    PCA_2D[:,1],
    s=40,
    c=alledaten.target
)

plt.show()