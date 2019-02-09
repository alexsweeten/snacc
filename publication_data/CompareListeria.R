require(ape)
require(phangorn)
require(treespace)

publication_trees<- readRDS("publication_trees.rds")

# reference tree obtained from https://github.com/johnlees/which_tree/blob/master/tree_compare.R
listeria_realtr <- midpoint(read.tree(paste(sep="/","benchmark_trees/RealTree_Listeria.nwk")))
listeria_realtr$edge.length <- listeria_realtr$edge.length*0.01 # correct for scaling introduced by ALF
listeria_samples <- sort(listeria_realtr$tip.label)

# Draw trees from distance matrices
temp = read.csv("listeria_distances/mash_listeria_distances.csv", sep=",")
listeria_mash.matrix <- as.matrix(temp, header=TRUE)
dimnames(listeria_mash.matrix) = list(listeria_samples, listeria_samples)

temp = read.csv("listeria_distances/listeria_snacc_lzma.csv", sep=",")
listeria_snacc_lzma.matrix <- as.matrix(temp)
dimnames(listeria_snacc_lzma.matrix) = list(listeria_samples, listeria_samples)

temp = read.csv("listeria_distances/listeria_snacc_lzma_reverse.csv", sep=",")
listeria_snacc_lzma_reverse.matrix <- as.matrix(temp)
dimnames(listeria_snacc_lzma_reverse.matrix) = list(listeria_samples, listeria_samples)

temp = read.csv("listeria_distances/listeria_snacc_lz4.csv", sep=",")
listeria_snacc_lz4.matrix <- as.matrix(temp)
dimnames(listeria_snacc_lz4.matrix) = list(listeria_samples, listeria_samples)

temp = read.csv("listeria_distances/listeria_snacc_bzip2.csv", sep=",")
listeria_snacc_bzip2.matrix <- as.matrix(temp)
dimnames(listeria_snacc_bzip2.matrix) = list(listeria_samples, listeria_samples)

temp = read.csv("listeria_distances/listeria_snacc_gzip.csv", sep=",")
listeria_snacc_gzip.matrix <- as.matrix(temp)
dimnames(listeria_snacc_gzip.matrix) = list(listeria_samples, listeria_samples)

temp = read.csv("listeria_distances/listeria_snacc_zlib.csv", sep=",")
listeria_snacc_zlib.matrix <- as.matrix(temp)
dimnames(listeria_snacc_zlib.matrix) = list(listeria_samples, listeria_samples)

listeria_andi_bionj <- publication_trees[["BIONJ + andi dist"]]

temp = read.csv("listeria_distances/poppunk_listeria_13.csv", sep=",")
listeria_poppunk13.matrix <- as.matrix(temp)
dimnames(listeria_poppunk13.matrix) = list(listeria_samples, listeria_samples)

temp = read.csv("listeria_distances/poppunk_listeria_17.csv", sep=",")
listeria_poppunk17.matrix <- as.matrix(temp)
dimnames(listeria_poppunk17.matrix) = list(listeria_samples, listeria_samples)

temp = read.csv("listeria_distances/poppunk_listeria_21.csv", sep=",")
listeria_poppunk21.matrix <- as.matrix(temp)
dimnames(listeria_poppunk21.matrix) = list(listeria_samples, listeria_samples)

temp = read.csv("listeria_distances/poppunk_listeria_25.csv", sep=",")
listeria_poppunk25.matrix <- as.matrix(temp)
dimnames(listeria_poppunk25.matrix) = list(listeria_samples, listeria_samples)

temp = read.csv("listeria_distances/poppunk_listeria_29.csv", sep=",")
listeria_poppunk29.matrix <- as.matrix(temp)
dimnames(listeria_poppunk29.matrix) = list(listeria_samples, listeria_samples)

listeria_mash_bionj <- midpoint(bionj(listeria_mash.matrix))
listeria_mash_upgma <- midpoint(upgma(listeria_mash.matrix))
listeria_snacc_lzma_bionj <- midpoint(bionj(listeria_snacc_lzma.matrix))
listeria_snacc_lzma_upgma <- midpoint(upgma(listeria_snacc_lzma.matrix))
listeria_snacc_lzma_reverse_bionj <- midpoint(bionj(listeria_snacc_lzma_reverse.matrix))
listeria_snacc_lzma_reverse_upgma <- midpoint(upgma(listeria_snacc_lzma_reverse.matrix))
listeria_snacc_lz4_bionj <- midpoint(bionj(listeria_snacc_lz4.matrix))
listeria_snacc_lz4_upgma <- midpoint(upgma(listeria_snacc_lz4.matrix))
listeria_snacc_bzip2_bionj <- midpoint(bionj(listeria_snacc_bzip2.matrix))
listeria_snacc_bzip2_upgma <- midpoint(upgma(listeria_snacc_bzip2.matrix))
listeria_snacc_gzip_bionj <- midpoint(bionj(listeria_snacc_gzip.matrix))
listeria_snacc_gzip_upgma <- midpoint(upgma(listeria_snacc_gzip.matrix))
listeria_snacc_zlib_bionj <- midpoint(bionj(listeria_snacc_zlib.matrix))
listeria_snacc_zlib_upgma <- midpoint(upgma(listeria_snacc_zlib.matrix))

listeria_poppunk13_bionj <- midpoint(bionj(listeria_poppunk13.matrix))
listeria_poppunk13_upgma <- midpoint(upgma(listeria_poppunk13.matrix))
listeria_poppunk17_bionj <- midpoint(bionj(listeria_poppunk17.matrix))
listeria_poppunk17_upgma <- midpoint(upgma(listeria_poppunk17.matrix))
listeria_poppunk21_bionj <- midpoint(bionj(listeria_poppunk21.matrix))
listeria_poppunk21_upgma <- midpoint(upgma(listeria_poppunk21.matrix))
listeria_poppunk25_bionj <- midpoint(bionj(listeria_poppunk25.matrix))
listeria_poppunk25_upgma <- midpoint(upgma(listeria_poppunk25.matrix))
listeria_poppunk29_bionj <- midpoint(bionj(listeria_poppunk29.matrix))
listeria_poppunk29_upgma <- midpoint(upgma(listeria_poppunk29.matrix))

#Create identity vectors
test_listeria <- diag(listeria_snacc_bzip2.matrix)
listeria_bzip2_identity <- mean(test_listeria)
test_listeria <- diag(listeria_snacc_gzip.matrix)
listeria_gzip_identity <- mean(test_listeria)
test_listeria <- diag(listeria_snacc_lz4.matrix)
listeria_lz4_indentity <- mean(test_listeria)
test_listeria <- diag(listeria_snacc_lzma.matrix)
listeria_lzma_identity <- mean(test_listeria)
test_listeria <- diag(listeria_snacc_zlib.matrix)
listeria_zlib_identity <- mean(test_listeria)

#Create correlation matrix

matrix_list_listeria <- list(mash = dist(listeria_mash.matrix),
                    bzip2 = dist(listeria_snacc_bzip2.matrix),
                    gzip = dist(listeria_snacc_gzip.matrix),
                    lz4 = dist(listeria_snacc_lz4.matrix),
                    lzma = dist(listeria_snacc_lzma.matrix),
                    zlib = dist(listeria_snacc_zlib.matrix))

temp <- c()
for(x in matrix_list_listeria){
  for(y in matrix_list_listeria){
    z <- mantel.rtest(x,y,nrepet=100)
    temp <- c(temp, z)
  }
}

correlation_lis.matrix <- matrix(
  temp,
  nrow=6,
  ncol=6
)
row.names(correlation_lis.matrix) <- names(matrix_list_listeria)
colnames(correlation_lis.matrix) <- names(matrix_list_listeria)

#plot correlation matrix
test <- corrplot(correlation_lis.matrix, type = "lower", order = "hclust", 
         tl.col = "black", tl.srt = 35, insig="p-value", sig.level = -1)

plot(test)

#test <- cophylo(listeria_realtr, listeria_snacc_lzma_bionj, fsize=0.7)
#plot(test, fsize=0.5)
