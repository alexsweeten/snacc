import numpy as np
import scipy.spatial
import scipy.cluster
import matplotlib.pyplot as plt
plt.switch_backend('agg')

from misc import read_dist

def metrify(D):
  # Converts a raw distance matrix to one which
  # is symmetric and has zero diagonals
  D_sym = 0.5*(D + D.T)
  np.fill_diagonal(D_sym, 0.)
  return D_sym

def hierarchical(D_sym):
  # Produce a condensed distance matrix
  D_cond = scipy.spatial.distance.squareform(D_sym)
  linkage = scipy.cluster.hierarchy.linkage(D_cond, method='average')
  return linkage
  
def plot_hierarchical(linkage, plt_file,n):
  plt.figure(figsize=(min(3+0.1*n,13),4))
  dn = scipy.cluster.hierarchy.dendrogram(linkage)
  plt.savefig(plt_file, bbox_inches="tight")
    
def main(csv_file, plt_file):
  # Does everything in order
  D = read_dist(csv_file)
  D_sym = metrify(D)
  n,_ = D_sym.shape
  linkage = hierarchical(D_sym)
  plot_hierarchical(linkage, plt_file,n)
  
def _test_dist2tree():
  #main("../test_dataset/distance_matrix_mysteryGenome1-8.csv", "../test_dataset/distance_matrix_mysteryGenome1-8_tree.png")
  main("../test_dataset/ecoli_mash.csv", "../test_dataset/ecoli_mash_tree.png")
if __name__=="__main__":
  _test_dist2tree()