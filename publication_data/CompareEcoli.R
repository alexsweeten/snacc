require(ape)
require(phangorn)
require(treespace)

ecoli_realtr <- midpoint(read.tree(paste(sep="/","benchmark_trees/RealTree_Ecoli.nwk")))
ecoli_realtr$edge.length <- ecoli_realtr$edge.length*0.01 # correct for scaling introduced by ALF
ecoli_samples <- sort(ecoli_realtr$tip.label)

temp = read.csv("ecoli_st131_distances/ecoli_snacc.csv", sep=",")
ecoli_snacc.matrix <- as.matrix(temp, header=TRUE)
dimnames(ecoli_snacc.matrix) = list(ecoli_samples, ecoli_samples)

temp = read.csv("ecoli_st131_distances/ecoli_snacc_reverse.csv", sep=",")
ecoli_snacc_reverse.matrix <- as.matrix(temp, header=TRUE)
dimnames(ecoli_snacc_reverse.matrix) = list(ecoli_samples, ecoli_samples)

temp = read.csv("ecoli_st131_distances/ecoli_mash.csv", sep=",")
ecoli_mash.matrix <- as.matrix(temp, header=TRUE)
dimnames(ecoli_mash.matrix) = list(ecoli_samples, ecoli_samples)

temp = read.csv("ecoli_st131_distances/ecoli_andi.csv", sep=",")
ecoli_andi.matrix <- as.matrix(temp, header=TRUE)
dimnames(ecoli_andi.matrix) = list(ecoli_samples, ecoli_samples)

temp = read.csv("ecoli_st131_distances/ecoli_poppunk_15.csv", sep=",")
ecoli_poppunk_15.matrix <- as.matrix(temp, header=TRUE)
dimnames(ecoli_poppunk_15.matrix) = list(ecoli_samples, ecoli_samples)

temp = read.csv("ecoli_st131_distances/ecoli_poppunk_19.csv", sep=",")
ecoli_poppunk_19.matrix <- as.matrix(temp, header=TRUE)
dimnames(ecoli_poppunk_19.matrix) = list(ecoli_samples, ecoli_samples)

temp = read.csv("ecoli_st131_distances/ecoli_poppunk_23.csv", sep=",")
ecoli_poppunk_23.matrix <- as.matrix(temp, header=TRUE)
dimnames(ecoli_poppunk_23.matrix) = list(ecoli_samples, ecoli_samples)

temp = read.csv("ecoli_st131_distances/ecoli_poppunk_27.csv", sep=",")
ecoli_poppunk_27.matrix <- as.matrix(temp, header=TRUE)
dimnames(ecoli_poppunk_27.matrix) = list(ecoli_samples, ecoli_samples)

ecoli_mash_bionj <- midpoint(bionj(ecoli_mash.matrix))
ecoli_mash_upgma <- midpoint(upgma(ecoli_mash.matrix))
ecoli_snacc_bionj <- midpoint(bionj(ecoli_snacc.matrix))
ecoli_snacc_upgma <- midpoint(upgma(ecoli_snacc.matrix))
ecoli_snacc_reverse_bionj <- midpoint(bionj(ecoli_snacc_reverse.matrix))
ecoli_snacc_reverse_upgma <- midpoint(upgma(ecoli_snacc_reverse.matrix))
ecoli_andi_bionj <- midpoint(bionj(ecoli_andi.matrix))
ecoli_andi_upgma <- midpoint(upgma(ecoli_andi.matrix))


ecoli_poppunk_15_bionj <- midpoint(bionj(ecoli_poppunk_15.matrix))
ecoli_poppunk_15_upgma <- midpoint(upgma(ecoli_poppunk_15.matrix))
ecoli_poppunk_19_bionj <- midpoint(bionj(ecoli_poppunk_19.matrix))
ecoli_poppunk_19_upgma <- midpoint(upgma(ecoli_poppunk_19.matrix))
ecoli_poppunk_23_bionj <- midpoint(bionj(ecoli_poppunk_23.matrix))
ecoli_poppunk_23_upgma <- midpoint(upgma(ecoli_poppunk_23.matrix))
ecoli_poppunk_27_bionj <- midpoint(bionj(ecoli_poppunk_27.matrix))
ecoli_poppunk_27_upgma <- midpoint(upgma(ecoli_poppunk_27.matrix))

