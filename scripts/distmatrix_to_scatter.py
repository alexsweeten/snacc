from sklearn import manifold
from sklearn import cluster
from sklearn import random_projection
from sklearn.decomposition import PCA
import umap
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from misc import read_dist
  
def reduce_dimension(D, projection='mds'):
  projections = {
         # For use with Distance Matrices:
         'mds' : manifold.MDS(2, dissimilarity="precomputed"),
         'tsne' : manifold.TSNE(2, metric="precomputed"),
         'umap' : umap.UMAP(n_components=2, metric='precomputed'),
         # For use with Affinity Matrices:
         'spectralembedding' : manifold.SpectralEmbedding(2, affinity="precomputed"),
         # These should not be used:
         'pca' : PCA(2),
         'gaussianrp' : random_projection.GaussianRandomProjection(2),
  }
  X = projections[projection].fit_transform(D)
  return X
  
def clustering(M, algorithm='affinity', **kwargs):
  algorithms = {
      # For use with affinity matrices. Note that smaller affinity means more separation.
      'affinity' : sklearn.cluster.AffinityPropagation(affinity='precomputed', **kwargs),
      'spectral' : sklearn.manifold.SpectralEmbedding(affinity='precomputed', **kwargs),
      'agglomerative' : sklearn.cluster.AgglomerativeClustering(affinity='precomputed', **kwargs),
      # For use with projections:
      'kmeans' : sklearn.cluster.KMeans(**kwargs),
      'dbscan' : sklearn.cluster.dbscan(**kwargs),
  }
  return algorithms[algorithm].fit_predict(M)
  
def f_inverse(D):
  return ((1./D) - 1.)


def distance_to_affinity(D, function=f_inverse):
  return function(D)

def _test_plot(csv_file, plt_file):
  projections = ['mds', 'tsne', 'umap', 'spectralembedding', 'pca', 'gaussianrp']
  N = len(projections)
  D = read_dist(csv_file)
  fig, ax = plt.subplots(N,1, figsize=(6,N*3))
  for i,p in enumerate(projections):
    X = reduce_dimension(D, p)
    ax[i].scatter(*X.T)
    ax[i].set_title(p)
  plt.tight_layout()
  plt.savefig(plt_file, bbox_inches="tight")
  
if __name__ == "__main__":
  #_test_plot("../test_dataset/distance_matrix_mysteryGenome1-8_skew.csv", "../test_dataset/distance_matrix_mysteryGenome1-8_skew_scatter.png")
  _test_plot("../test_dataset/sample100.csv", "../test_dataset/sample100_skew_scatter.png")
