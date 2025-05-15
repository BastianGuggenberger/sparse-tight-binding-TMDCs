#Wir wollen das 4-D Iris-Dataset 3-D machen
from sklearn.datasets import load_iris
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

daten = load_iris(as_frame=True)
#print(daten.keys())
#print(daten['data'])
#print(daten['target'])
#print('FRAME:::')
#print(daten['frame'])



# Rename classes using the iris target names
#daten.frame["target"] = daten.target_names[daten.target]
#print('FRAME:::')
#print(daten['frame'])
#plot = sns.pairplot(daten.frame, hue="target")
#plt.show()



# unused but required import for doing 3d projections with matplotlib < 3.2
#import mpl_toolkits.mplot3d  # noqa: F401


#fig ist das Fenster
#ax ist der Plot der auf dem Fenster gezeigt wir
fig = plt.figure(1, figsize=(8, 6))
ax = fig.add_subplot(111)


print(daten.data)
X_reduced = PCA(n_components=2).fit_transform(daten.data)
print(X_reduced)
print(type(X_reduced))
scatter = ax.scatter(
    X_reduced[:, 0],
    X_reduced[:, 1],
    c=daten.target,
    s=40,
)
print(X_reduced[0,0])


ax.set(
    title="First three PCA dimensions",
    xlabel="1st Eigenvector",
    ylabel="2nd Eigenvector",
)
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])

# Add a legend
legend1 = ax.legend(
    scatter.legend_elements()[0],
    daten.target_names.tolist(),
    loc="upper right",
    title="Classes",
)
ax.add_artist(legend1)

plt.show()
