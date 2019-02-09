require(ape)
require(phangorn)
require(treespace)

#Maximum clade credibility tree obtained from the R package "treespace." Jombart et. al, 2018.
dengue_realtr <- DengueBEASTMCC
dengue_samples <- sort(dengue_realtr$tip.label)

#Import distance matrices from file
temp = read.csv("dengue_distances/mash_dengue.csv", sep=",")
mash_dengue.matrix <- as.matrix(temp, header=TRUE)
dimnames(mash_dengue.matrix) = list(dengue_samples, dengue_samples)

temp = read.csv("dengue_distances/snacc_bzip2_dengue.csv", sep=",")
snacc_bzip2_dengue.matrix <- as.matrix(temp, header=TRUE)
dimnames(snacc_bzip2_dengue.matrix) = list(dengue_samples, dengue_samples)

temp = read.csv("dengue_distances/snacc_gzip_dengue.csv", sep=",")
snacc_gzip_dengue.matrix <- as.matrix(temp, header=TRUE)
dimnames(snacc_gzip_dengue.matrix) = list(dengue_samples, dengue_samples)

temp = read.csv("dengue_distances/snacc_lz4_dengue.csv", sep=",")
snacc_lz4_dengue.matrix <- as.matrix(temp, header=TRUE)
dimnames(snacc_lz4_dengue.matrix) = list(dengue_samples, dengue_samples)

temp = read.csv("dengue_distances/snacc_lzma_dengue.csv", sep=",")
snacc_lzma_dengue.matrix <- as.matrix(temp, header=TRUE)
dimnames(snacc_lzma_dengue.matrix) = list(dengue_samples, dengue_samples)

temp = read.csv("dengue_distances/snacc_zlib_dengue.csv", sep=",")
snacc_zlib_dengue.matrix <- as.matrix(temp, header=TRUE)
dimnames(snacc_zlib_dengue.matrix) = list(dengue_samples, dengue_samples)

#Create list of matrices
matrix_list <- list(mash = mash_dengue.matrix,
                    bzip2 = snacc_bzip2_dengue.matrix,
                    gzip = snacc_gzip_dengue.matrix,
                    lz4 = snacc_lz4_dengue.matrix,
                    lzma = snacc_lzma_dengue.matrix,
                    zlib = snacc_zlib_dengue.matrix)

#Create trees from distance matrices
mash_dengue_tr <- root(bionj(mash_dengue.matrix), resolve.root=TRUE, outgroup="D4Thai63")
snacc_bzip2_dengue_tr <- root(bionj(snacc_bzip2_dengue.matrix), resolve.root=TRUE, outgroup="D4Thai63")
snacc_gzip_dengue_tr <- root(bionj(snacc_gzip_dengue.matrix), resolve.root=TRUE, outgroup="D4Thai63")
snacc_lz4_dengue_tr <- root(bionj(snacc_lz4_dengue.matrix), resolve.root=TRUE, outgroup="D4Thai63")
snacc_lzma_dengue_tr <- root(bionj(snacc_lzma_dengue.matrix), resolve.root=TRUE, outgroup="D4Thai63")
snacc_zlib_dengue_tr <- root(bionj(snacc_zlib_dengue.matrix), resolve.root=TRUE, outgroup="D4Thai63")

#Create correlation matrix
temp <- c()
for(x in matrix_list){
  for(y in matrix_list){
    z <- cor(c(x),c(y))
    temp <- c(temp, z)
  }
}

correlation.matrix <- matrix(
  temp,
  nrow=6,
  ncol=6
)
row.names(correlation.matrix) <- names(matrix_list)
colnames(correlation.matrix) <- names(matrix_list)

#plot correlation matrix
corrplot(correlation.matrix, type = "lower", order = "hclust", 
         tl.col = "black", tl.srt = 35, method = "number")

#Create list of trees
tree_list <- list( reference = dengue_realtr,
                   mash = mash_dengue_tr,
                   bzip2 = snacc_bzip2_dengue_tr,
                   gzip = snacc_gzip_dengue_tr,
                   lz4 = snacc_lz4_dengue_tr,
                   lzma = snacc_lzma_dengue_tr,
                   zlib = snacc_zlib_dengue_tr)

#Create matrix of Kendall Colijn (KC) distances
temp <- c()
for(x in tree_list){
  for(y in tree_list){
    z <- treeDist(x,y)
    temp <- c(temp, z)
  }
}

KC.matrix <- matrix(
  temp,
  nrow=7,
  ncol=7
)

row.names(KC.matrix) <- names(tree_list)
colnames(KC.matrix) <- names(tree_list)

#Create matrix of Robinson Fould (RF) distances
temp <- c()
for(x in tree_list){
  for(y in tree_list){
    z <- RF.dist(x,y)
    temp <- c(temp, z)
  }
}

RF.matrix <- matrix(
  temp,
  nrow=7,
  ncol=7
)

row.names(RF.matrix) <- names(tree_list)
colnames(RF.matrix) <- names(tree_list)

#Create matrix of Subtree Branch & Repruning (SPR) distances
temp <- c()
for(x in tree_list){
  for(y in tree_list){
    z <- SPR.dist(x,y)
    temp <- c(temp, z)
  }
}

SPR.matrix <- matrix(
  temp,
  nrow=7,
  ncol=7
)

row.names(SPR.matrix) <- names(tree_list)
colnames(SPR.matrix) <- names(tree_list)

#Create vector of identity, Symmetry, and triangle inequality properties.
test_dengue <- diag(snacc_lzma_dengue.matrix)
dengue_lzma_identity <- mean(test_dengue)
test_dengue <- diag(snacc_lz4_dengue.matrix)
dengue_lz4_identity <- mean(test_dengue)
test_dengue <- diag(snacc_bzip2_dengue.matrix)
dengue_bzip2_identity <- mean(test_dengue)
test_dengue <- diag(snacc_gzip_dengue.matrix)
dengue_gzip_identity <- mean(test_dengue)
test_dengue <- diag(snacc_zlib_dengue.matrix)
dengue_zlib_identity <- mean(test_dengue)
test_dengue <- diag(mash_dengue.matrix)
mash_identity <- mean(test_dengue)
