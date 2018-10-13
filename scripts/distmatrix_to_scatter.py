from sklearn import manifold
from sklearn import cluster
from sklearn import random_projection
from sklearn.decomposition import PCA
import umap
import numpy as np
import matplotlib.pyplot as plt
from misc import read_dist
plt.switch_backend('agg')

color_palette = plt.cm.Set1


def reduce_dimension(D, projection='mds'):
    projections = {
                 # For use with Distance Matrices:
                 'mds': manifold.MDS(2, dissimilarity="precomputed"),
                 'tsne': manifold.TSNE(2, metric="precomputed"),
                 'umap': umap.UMAP(n_components=2, metric='precomputed'),
                 # For use with Affinity Matrices:
                 'spectralembedding': manifold.SpectralEmbedding(2, affinity="precomputed"),
                 # These should not be used:
                 'pca': PCA(2),
                 'gaussianrp': random_projection.GaussianRandomProjection(2)
    }
    X = projections[projection].fit_transform(D)
    return X


def clustering(M, algorithm='agglomerative', **kwargs):
    algorithms = {
            # For use with affinity matrices. Note that smaller affinity means more separation.
            'affinity': cluster.AffinityPropagation(affinity='precomputed', **kwargs),
            'spectral': cluster.SpectralClustering(affinity='precomputed', **kwargs),
            'agglomerative': cluster.AgglomerativeClustering(affinity='precomputed', linkage='average', **kwargs),
            # For use with projections:
            'kmeans': cluster.KMeans(M, **kwargs),
            'dbscan': cluster.dbscan(M, **kwargs)
    }
    return algorithms[algorithm].fit_predict(M)


def f_inverse(D):
    E = ((1./D) - 1.)
    E[np.logical_not(np.isfinite(E))] = 1e12
    np.fill_diagonal(E, 0.)
    return E


def distance_to_affinity(D, function=f_inverse):
    return function(D)


def plot_labels(ax, X, labels):
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    cluster_colors = [color_palette(x) if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in labels]
    ax.scatter(*X.T, c=cluster_colors)


def _test_all(csv_file, img_file):
    projections_dist = ['mds', 'tsne', 'umap']
    projections_affinity = ['spectralembedding']
    projections = projections_dist.copy()
    projections.extend(projections_affinity)
    algorithms = ['affinity', 'spectral', 'agglomerative']
    D = read_dist(csv_file)
    E = distance_to_affinity(D)
    Xs = [reduce_dimension(D, projection=p) for p in projections_dist]
    Xs_a = [reduce_dimension(E, projection=p) for p in projections_affinity]
    Xs.extend(Xs_a)
    Ls = [clustering(E, algorithm=a) for a in algorithms]
    M = len(algorithms)
    N = len(projections)
    fig, ax = plt.subplots(N, M, figsize=(M*3, N*3))
    for i, x in enumerate(Xs):
        for j, labels in enumerate(Ls):
            plot_labels(ax[i, j], x, labels)
            ax[i, j].set_title('{} {}'.format(projections[i], algorithms[j]))
    plt.tight_layout()
    plt.savefig(img_file, bbox_inches="tight")


if __name__ == "__main__":
    img_file = "../test_dataset/ecoli_mash_scatter.png"
    csv_file = "../test_dataset/ecoli_mash.csv"
    _test_all(csv_file, img_file)
