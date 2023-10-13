import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from typing import Union

plt.set_loglevel("error")

MARKERS = list(Line2D.markers.keys())


def get_best_ncluster(coords):
    """
    Obtain the best number of cluster based on the maximization of the silhouette

    :param coords: array of the coordinates of each atom
    :type coords: 2D-array

    :return: best number of clusters
    :rtype: int
    """

    k_range = range(10, 30)
    silhouette_scores = []
    for k in k_range:
        kmeans = KMeans(n_clusters=k, n_init=10)
        labels = kmeans.fit_predict(coords)
        score = silhouette_score(coords, labels)
        silhouette_scores.append(score)
    return k_range[np.argmax(silhouette_scores)]


def calc_pca(confs: list, ncluster: Union[int, None] = None) -> tuple:
    """
    Function that execute the actual PCA analysis.
    It wants to understand how conformations differ from each other based on their overall Cartesian coordinates

    :param confs: whole list of the confomers
    :type confs: list
    :param ncluster: number of cluster to form using the KMean analysis
    :type ncluster: int
    :return: PCA transformation, Clustered coordinates
    :rtype: tuple
    """

    # fetch all geometries and reshaping them to create the correct 2D matrix
    data = np.array([conf.last_geometry for conf in confs])
    colors = [conf.color for conf in confs]
    numbers = [conf.number for conf in confs]
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

    # Cluster the data the first time
    if not confs[0].cluster:
        kmeans = KMeans(n_clusters=n_c, n_init=10)
        clusters = kmeans.fit_predict(pca_scores)
        for idx, conf in enumerate(confs):
            conf.cluster = clusters[idx]
    else:
        clusters = [conf.cluster for conf in confs]

    return pca_scores, clusters, colors, numbers


def obtain_markers_from_cluster(cluster: int):
    """
    Obtain a different marker from the marker library for different conformers

    :param cluster: the cluster number
    :return: marker
    :rtype: matplotlib.lines
    """
    return MARKERS[cluster]


def save_PCA_snapshot(
    fname: str,
    title: str,
    pca_scores: np.ndarray,
    clusters: np.ndarray,
    colors: list,
    numbers: list,
):
    """
    Graph and save the image of the PCA analysis

    :param fname: filename to save the graphs
    :type fname: str
    :param title: title of the graph
    :type title: str
    :param pca_scores: PCA transformation
    :type pca_scores: np.array
    :param clusters: Clustered coordinates
    :type clusters: np.array
    :rtype: None
    """

    fig, ax = plt.subplots()

    for x, y, m, c, n in zip(
        pca_scores[:, 0],
        pca_scores[:, 1],
        np.array(list(map(obtain_markers_from_cluster, clusters))),
        colors,
        numbers,
    ):
        ax.scatter(x, y, c=c, marker=m, label=f"CONF {n}")
    ax.set_xlabel("Principal Component 1")
    ax.set_ylabel("Principal Component 2")
    ax.set_title(title)

    y = ax.get_ylim()
    x = ax.get_xlim()

    ax.vlines(0, y[0], y[1], "#353535", "--", alpha=0.2)
    ax.hlines(0, x[0], x[1], "#353535", "--", alpha=0.2)

    ax.set_xlim(x)
    ax.set_ylim(y)
    plt.legend(
        loc="upper left",
        bbox_to_anchor=(1.05, 1.0),
        fancybox=True,
        shadow=True,
        ncol=2,
        title="Conformers",
    )
    plt.tight_layout()
    plt.savefig(fname, dpi=300)
    return None


def perform_PCA(confs: list, ncluster: int, fname: str, title: str, log) -> None:
    """
    Perform a PCA analysis

    :param confs:  list of all the active conformers
    :type confs: list
    :param ncluster:  number of cluster to group the ensemble
    :type ncluster: int
    :param fname:  filename for the graph
    :type fname: str
    :param title:  title of the graph
    :type title: str
    :param log:  logger instance
    :type log: logging
    :rtype: None
    """
    log.info("Starting PCA analysis")
    nc = ncluster if len(confs) > ncluster else len(confs) - 1
    if nc <= 2:
        return None
    pca_scores, clusters, colors, numbers = calc_pca(confs, ncluster=nc)
    save_PCA_snapshot(fname, title, pca_scores, clusters, colors, numbers)

    return None


if __name__ == "__main__":  # pragma: no cover:
    from ioFile import read_ensemble
    import mock

    # Load the XYZ file
    xyz_file = read_ensemble("files/ensemble.xyz", 0, 1, mock.MagicMock())
    perform_PCA(xyz_file, 5, "files/test.png", "Test", mock.MagicMock())
