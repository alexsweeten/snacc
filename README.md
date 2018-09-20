# bioncd-hackseq
**BIO**logical implementation of the **N**ormalized **C**ompression **D**istance (name tentative to change). This repository was initated for the [hackseq18](https://www.hackseq.com) project titled "Alignment-free Pathogen Genomics", with development scheduled for October 12th - 14th, 2018.

# Initial Description
Sequence alignment has long been the de-facto method for determing similarity between two genomes. Algorithms which implement sequence alignment are abundant and widely used; the Basic Local Alignment Search Tool algorithm (BLAST) is one of the most frequently cited scientific papers. Recent advances in sequencing technology and bioinformatics toolkits have resulted in an explosive growth of genomic data available to analyze. iven the explosion in the amount of genomic data over hte past decade, genomic information has become as well as data stored in online databases requires alternative methods.

The normalized compression distance is a paramater-free distance metric, which can in be applied towards genetic sequences.
Intuition states that the more similar two sequences are, the It is a real-world approximation of Kolmogorov Complexity 

# Project Goals 
Our goal is to develop BioNCD, a pipeline that takes sequences as input, and outputs a sequence similairity network with NCD used as a distance metric. This pipeline can be broken down into 3 stages:

1) Compression
* Input: Set of sequences
* Output: Distance matrix of NCD values

2) Clustering
* Input: Distance matrix of NCD values
* Output: Table of cluster assignments

3) Visualization
* Input: Table of cluster assignments & distance matrix of NCD values
* Output: Sequence similarity network

The compression stage of this pipeline is where we will focus the bulk of our development towards. Finally, we would like to package this pipeline into a conda package to make it easy for users to run this pipeline and not have any dependency issues. 

# Hackathon pre-requisities & structure

I have initialized this repository with a number of issues for participants to work on. As work gets completed, I will add more issues as they become available.

For communication, we plan to use a Slack channel. 
If you are uncomfortable, don't fret! We 
Knowledge of compression algorithms/information theory is not required, but would be nice to know. 

# Dependencies
Dependencies are not finalized and will be updated here. Most likely we will require:
- [`python3`](https://python.org)
- [`numpy`](https://numpy.org)
- [`scipy`](https://scipy.org)
- [`pandas`](https://pandas.pydata.org)
- [`matplotlib`](https://matplotlib.org)
- [`sklearn`](http://scikit-learn.org/stable/)

