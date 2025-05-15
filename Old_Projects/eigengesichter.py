#https://scikit-learn.org/stable/auto_examples/applications/plot_face_recognition.html#sphx-glr-auto-examples-applications-plot-face-recognition-py

from time import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import loguniform

from sklearn.datasets import fetch_lfw_people
from sklearn.decomposition import PCA
from sklearn.metrics import ConfusionMatrixDisplay, classification_report
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


lfw_people = fetch_lfw_people(min_faces_per_person=70, resize=0.4)

images = lfw_people.images
data = lfw_people.data
targets = lfw_people.target
targetnames = lfw_people.target_names

print(type(lfw_people))
print(type(images),type(data),type(targets),type(targetnames))
print(images.shape,data.shape,targets.shape,targetnames.shape)

print(images)
print(data)
print(targets)
print(targetnames)
print(data[0])


#reduziert_10d = PCA(n_components=10, svd_solver='randomized').fit(data)
#eigenvektoren_10d = reduziert_10d.components_
#reduziert_10d_daten = PCA(n_components=10, svd_solver='randomized').fit_transform(data)

#print(eigenvektoren_10d)

def ruecktransformation(data,dimensions):
    neuekords = PCA(n_components=dimensions,svd_solver='randomized').fit_transform(data)
    neuebasis = PCA(n_components=dimensions,svd_solver='randomized').fit(data).components_
    #print(neuekords.shape)
    #print(neuebasis.shape)
    return np.matmul(neuekords,neuebasis)

def plot_gallery(exakt, h, w, n_row=5, n_col=4, steps=30):
    plt.figure(figsize=(1.8 * n_col, 2.4 * n_row))
    plt.subplots_adjust(bottom=0, left=0.01, right=0.99, top=0.90, hspace=0.35)
    for i in range(0,n_row):
        print(i)
        for j in range(n_col):
            plt.subplot(n_row, n_col, i*n_col + j + 1)
            approximiert = ruecktransformation(exakt,1+i*steps)
            plt.imshow(approximiert[j].reshape((h, w)), cmap=plt.cm.gray)
            titel = targetnames[targets[j]] + " \n " + str(1+i*steps) + " dimensions"
            plt.title( titel, size=12)
            plt.xticks(())
            plt.yticks(())
    for j in range(n_col):
        plt.subplot(n_row, n_col, i*n_col + j + 1)
        plt.imshow(exakt[j].reshape((h, w)), cmap=plt.cm.gray)
        titel = targetnames[targets[j]] + "\n exakt"
        plt.title(titel, size=12)
        plt.xticks(())
        plt.yticks(())



#plot_gallery(data, images.shape[1], images.shape[2])
plt.show()
