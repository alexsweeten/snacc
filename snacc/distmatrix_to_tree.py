import scipy.spatial
import scipy.cluster
import matplotlib.pyplot as plt
from misc import read_dist_values_names
from misc import metrify
plt.switch_backend('agg')


def hierarchical(D_sym):
    # Produce a condensed distance matrix
    D_cond = scipy.spatial.distance.squareform(D_sym)
    # Calculate the tree
    linkage = scipy.cluster.hierarchy.linkage(D_cond, method='average')
    return linkage


def plot_hierarchical(linkage, plt_file, n, labels=None):
    plt.figure(figsize=(min(3+0.1*n, 13), 4))
    dn = scipy.cluster.hierarchy.dendrogram(linkage, labels=labels)
    plt.savefig(plt_file, bbox_inches="tight", dpi=300)


def get_newick(node, newick, parentdist, leaf_names):
    # Credit: user jfn on stackoverflow
    # stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id],
                              parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), newick, node.dist, leaf_names)
        newick = get_newick(node.get_right(), ",%s" % (newick),
                           node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick


def write_newick(linkage, leaf_names, nwk_file):
    tree = scipy.cluster.hierarchy.to_tree(linkage, False)
    nwk_string = get_newick(tree, "", tree.dist, leaf_names)
    with open(nwk_file, "w+") as f:
        f.write(nwk_string)


def main(csv_file, plt_file, nwk_file):
    # Does everything in order
    leaf_names, D = read_dist_values_names(csv_file)
    D_sym = metrify(D)
    n, _ = D_sym.shape
    linkage = hierarchical(D_sym)
    plot_hierarchical(linkage, plt_file, n, labels=leaf_names)
    write_newick(linkage, leaf_names, nwk_file)


def _test_dist2tree():
    main("../test_dataset/streptococcus_lzma_dist_matrix.csv",
         "../test_dataset/streptococcus_lzma_tree.png",
         "../test_dataset/streptococcus_lzma_tree.nwk")


if __name__ == "__main__":
    _test_dist2tree()
