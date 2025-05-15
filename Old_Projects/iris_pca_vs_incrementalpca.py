#Vergleich von gewöhnlicher PCA mit Incremental PCA

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from sklearn.decomposition import PCA, IncrementalPCA


iris = load_iris()
rohdaten = iris.data
targets = iris.target


daten_ipca = IncrementalPCA(n_components=2, batch_size=10).fit_transform(rohdaten)
daten_pca = PCA(n_components=2).fit_transform(rohdaten)

print(type(daten_ipca))
print(daten_ipca.shape)
print(type(daten_pca))
print(daten_pca.shape)

ax = plt.axes(projection = "3d")

for datensatz, hoehe in [(daten_ipca, 0.0),(daten_pca,10.0)]:
    ax.scatter(
        daten_ipca[:, 0],
        daten_ipca[:, 1],
        hoehe
    )
    


"""
colors = ["navy", "turquoise", "darkorange"]
for X_transformed, title in [(daten_ipca, "Incremental PCA"), (daten_pca, "PCA")]:
    plt.figure(figsize=(8, 8))
    for color, i, target_name in zip(colors, [0, 1, 2], iris.target_names):
        plt.scatter(
            X_transformed[targets == i, 0],
            X_transformed[targets == i, 1],
            color=color,
            lw=2,
            label=target_name,
        )

    plt.title(title + " of iris dataset")
    plt.legend(loc="best", shadow=False, scatterpoints=1)
    plt.axis([-4, 4, -1.5, 1.5])
"""


plt.show()