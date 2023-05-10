import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from typing import Union


def get_best_ncluster(coords):
    """
    Obtain the best number of cluster based on the maximization of the silhouette 
    
    coords | 2D-array : array of the coordinates of each atom
    return | int : best number of clusters
    """

    k_range = range(10, 30)
    silhouette_scores = []
    for k in k_range:
        kmeans = KMeans(n_clusters=k, n_init=10)
        labels = kmeans.fit_predict(coords)
        score = silhouette_score(coords, labels)
        silhouette_scores.append(score)
    return k_range[np.argmax(silhouette_scores)]

def calc_pca(confs: list, ncluster : Union[int,None] = None) -> tuple:
    """
    Function that execute the actual PCA analysis. It wants to understand how conformations differ from each other based on their overall Cartesian coordinates

    confs | list : whole list of the confomers
    ncluster | int : number of cluster to form using the KMean analysis
    return tuple(np.array, np.array) : PCA transformation, Clustered coordinates
    """
    # fetch all geometries and reshaping them to create the correct 2D matrix
    data = np.array([atom.last_geometry for atom in confs])
    coords_2d = np.reshape(data, (data.shape[0], data.shape[1] * data.shape[2]))

    # normalize it to have mean=0 and variance=1
    coords_2d = (coords_2d - np.mean(coords_2d, axis=0)) / np.std(coords_2d, axis=0)

    # perform PCA analysis with number of components as minimum between number of 
    # n. confs and whole geom
    pca = PCA(n_components=min(coords_2d.shape[0], coords_2d.shape[1]))
    pca.fit(coords_2d)
    pca_scores = pca.transform(coords_2d)

    # getting the best number of clusters
    if not ncluster:
        n_c = get_best_ncluster(coords_2d)
    else:
        n_c = ncluster

    # Cluster the data
    kmeans = KMeans(n_clusters=n_c, n_init=10)
    clusters = kmeans.fit_predict(pca_scores)
    return pca_scores, clusters

def save_PCA_snapshot(fname : str, title : str, pca_scores : np.ndarray, clusters : np.ndarray) -> None: 
    """
    Graph and save the image of the PCA analysis
    
    fname | str : filename to save the graphs
    title | str : title of the graph
    pca_scores | np.array : PCA transformation
    clusters | np.array : Clustered coordinates
    return None
    """
    
    fig, ax = plt.subplots()
    ax.scatter(pca_scores[:,0], pca_scores[:,1], c=clusters, marker='o')
    ax.set_xlabel('Principal Component 1')
    ax.set_ylabel('Principal Component 2')
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(fname, dpi=300)
    return None


def perform_PCA(confs : list, ncluster : int, fname : str, title : str, log) -> None: 
    log.info('Starting PCA analysis')
    nc = ncluster if len(confs) > ncluster else len(confs)-1
    pca_scores, clusters = calc_pca(confs, ncluster=nc)
    save_PCA_snapshot(fname, title, pca_scores, clusters)

    return None

if __name__ == '__main__': 

    from ioFile import read_ensemble
    import mock

    # Load the XYZ file
    xyz_file = read_ensemble('files/ensemble.xyz', 0, 1, mock.MagicMock())
    perform_PCA(xyz_file, 5, 'files/test.png', 'Test', mock.MagicMock())