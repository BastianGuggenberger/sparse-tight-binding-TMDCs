#PCA der Gesichter auf 2 Dimensionen

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.datasets import fetch_lfw_people
from sklearn.decomposition import PCA

people = fetch_lfw_people()
images = people.images
data = people.data
targets = people.target
targetnames = people.target_names



##################
###TRAINING#######
##################

testset_list = []
for i in range(10,data.shape[0]):
    testset_list.append(data[i])
testset = np.array(testset_list)

images_2d = PCA(n_components=2, svd_solver='randomized').fit_transform(testset_list)
#print(data.shape)
#print(images_2d.shape)
pca = PCA(n_components=min(testset.shape),svd_solver='randomized').fit(testset_list)
eigengesichter = pca.components_
#print(eigengesichter.shape)

height = images.shape[1]
width = images.shape[2]
transformed_test_data = pca.transform(data)

data_backtransform = np.matmul(transformed_test_data,eigengesichter)


#13
n = 6
plt.figure(figsize=(1.8 * n, 2.4 * 1))
for i in range(0,n):
    plt.subplot(1, n, i+1)
    plt.imshow(images[i].reshape((height, width)), cmap=plt.cm.gray)
plt.savefig(fname="Gesichter_echt.jpg")