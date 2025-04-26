#Das Programm führt eine PCA Analyse durch, streicht n_reduzieren Eigenvektoren weg und transformiert dann die Daten zurück
#Danach werden die neuen Daten durch die Frobeniusnorm mit den alten Daten verglichen.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from sklearn.decomposition import PCA
import math


#FUNKTION frobeniusnorm: Berechnet die Frobeniusnorm einer Matrix
def frobeniusnorm(matrix):
    erg = 0
    for i in range(0, matrix.shape[0]-1):
        for j in range(0, matrix.shape[1]-1):
            eintrag= matrix[i][j]
            erg += eintrag * eintrag
    erg = math.sqrt(erg)
    return erg


#FUNKTION reduktionsfehler: Vergleich der Rücktransformierten nach Reduktion um n Komponenten
def reduktionsfehler(n_reduzieren, printvergleichsmatrix, printevs):

    # 1. PCA DURCHFÜHREN
    Rohdaten = load_iris().data
    Transformiert = PCA(n_components = 4-n_reduzieren).fit(Rohdaten) #Objekt aus dem wir später EVs extrahieren können
    Transformiert_Daten = PCA(n_components = 4-n_reduzieren).fit_transform(Rohdaten) #Transformierte Daten (Neue Koordinatenvektoren in Zeilen)
    evs = Transformiert.components_ #Neue Basis (Neue Basisvektoren in Zeilen)
    if (printevs == True):
        print('Eigenvektoren(n_reduzieren =', n_reduzieren, '): \n', evs)

    # 2. RÜCKTRANSFORMATION
    Ruecktransformiert = (np.matmul(evs.transpose(), Transformiert_Daten.transpose())).transpose()
    Ruecktransformiert += Transformiert.mean_
    #print(Ruecktransformiert)

    # 3. VERGLEICH
    vergleichsmatrix = np.empty(Ruecktransformiert.shape)
    for i in range(0, vergleichsmatrix.shape[0]-1):
        for j in range(0, vergleichsmatrix.shape[1]-1):
            vergleichsmatrix[i][j]=abs(Ruecktransformiert[i][j]-Rohdaten[i][j])
    if (printvergleichsmatrix == True):
        print(vergleichsmatrix)

    return frobeniusnorm(vergleichsmatrix)




# MAIN:
for i in range (5):
    print('Frobeniusnorm(n_reduzieren =',i,') =', reduktionsfehler(i,False,True), '\n')