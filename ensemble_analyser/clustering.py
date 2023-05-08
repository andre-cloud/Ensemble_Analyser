

import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def get_best_ncluster(coords):
    """
    Obtain the best number of cluster based on the maximization of the silhouette 
    
    return | int : best number of clusters"""
    k_range = range(10, 30)
    silhouette_scores = []
    for k in k_range:
        kmeans = KMeans(n_clusters=k, n_init=10)
        labels = kmeans.fit_predict(coords)
        score = silhouette_score(coords, labels)
        silhouette_scores.append(score)
    return k_range[np.argmax(silhouette_scores)]

def calc_pca(coords, ncluster = None):
    n_atoms = coords.shape[1]
    coords_2d = np.reshape(coords, (coords.shape[0], n_atoms * 3))


    pca = PCA(n_components=10)
    pca.fit(coords_2d)
    pca_scores = pca.transform(coords_2d)

    # getting the best number of clusters
    if not ncluster:
        optimal_k = get_best_ncluster(coords_2d)
    else:
        optimal_k = ncluster

    # Cluster the data
    kmeans = KMeans(n_clusters=optimal_k, n_init=10)
    clusters = kmeans.fit_predict(pca_scores)
    return pca_scores, clusters


if __name__ == '__main__': 
    
    from ioFile import read_ensemble

    # Load the XYZ file
    xyz_file = read_ensemble('files/ensemble.xyz', 0, 1, None)
    coords = np.array([atom.last_geometry for atom in xyz_file])
    
    pca_scores, clusters = calc_pca(coords, ncluster=5)

    import matplotlib.pyplot as plt 

    plt.style.context('bmh')
    plt.scatter(pca_scores[:,0], pca_scores[:,1], c=clusters)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.show()

    plt.show()