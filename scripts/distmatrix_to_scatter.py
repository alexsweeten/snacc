from sklearn import manifold
from sklearn import random_projection
from sklearn.decomposition import PCA
from sklearn import 
import umap
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from misc import read_dist
  
def reduce_dimension(D, projection='mds'):
  projections = {'mds' : manifold.MDS(2, dissimilarity="precomputed"),
                 'tsne' : manifold.TSNE(2, metric="precomputed"),
                 'gaussianrp' : random_projection.GaussianRandomProjection(2),
                 'spectralembedding' : manifold.SpectralEmbedding(2),
                 'pca' : PCA(2),
                 'umap' : umap.UMAP(n_components=2, metric='precomputed')
  }
  
  X = projections[projection].fit_transform(D)
  return X
  
def test_plot(csv_file, plt_file):
  projections = ['mds', 'tsne', 'gaussianrp', 'spectralembedding', 'pca', 'umap']
  N = len(projections)
  D = read_dist(csv_file)
  fig, ax = plt.subplots(N,1, figsize=(6,N*3))
  for i,p in enumerate(projections):
    X = reduce_dimension(D, p)
    ax[i].scatter(*X.T)
    ax[i].set_title(p)
  plt.tight_layout()
  plt.savefig(plt_file, bbox_inches="tight")
  
  
  
def clustering(csv_file, algorithm='kmeans', projection='mds'):
  D =  read_dist(csv_file)
  X = reduce_dimension(D, projection)
  
  if algorithm == 'kmeans':
    kmeans = cluster.KMeans(init='k-means++', n_clusters=3)
    y = kmeans.fit(X)
  
  if algorithm == 'agglomerative':
    agglomerative = cluster.AgglomerativeClustering(n_clusters=3)
    y = agglomerative.fit(X)
    
  return y
    
  
  
if __name__ == "__main__":
  #test_plot("../test_dataset/distance_matrix_mysteryGenome1-8_skew.csv", "../test_dataset/distance_matrix_mysteryGenome1-8_skew_scatter.png")
  test_plot("../test_dataset/sample100.csv", "../test_dataset/sample100_skew_scatter.png")
