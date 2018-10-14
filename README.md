# bioncd-hackseq
**BIO**logical implementation of the **N**ormalized **C**ompression **D**istance (name tentative to change). This repository was initated for the [hackseq18](https://www.hackseq.com) project titled "Alignment-free Pathogen Genomics", with development scheduled for October 12th - 14th, 2018.

## Initial Description
Sequence alignment has long been the de-facto method for determing similarity between two genomes. Algorithms which implement sequence alignment are abundant and widely used; the paper describing BLAST (Basic Local Alignment Search Tool algorithm) has over 70,000 citations, making it one of the most frequently cited tools in biology. However, in cases of low sequence homology, horizontal gene transfer, or lack of apriori information, as is common when dealing with pathogenic bacteria, alignment-based methods pose significant problems. New methods are required to analyze this data, and I hope to introduce the Normalized Compression Distance (NCD) as one such method.

NCD is a parameter-free metric which utilizes compression algorithms in order to estimate similarity. These algorithms find and take advantage of redundancies in a given input signal in order to reduce the size required to encode that signal. Intuition says that the more similar two strings are, the more redundancies they will have, which will result in a better compression score. NCD has been used in many different ways, such as [classifying musical genres of mp3 files](https://homepages.cwi.nl/~paulv/papers/music.pdf), [scanning Android files for viruses](https://link.springer.com/article/10.1007/s11416-015-0260-0) and [natural language processing](http://www.aclweb.org/anthology/P10-2015). However, there are currently no existing tools which apply NCD towards biological datasets.

## Project Goals
Our goal is to develop BioNCD, a pipeline that implements the normalized compression distance specifically for biological data as input. BioNCD will take sequences as input, and output a sequence similairity network with NCD used as the similarity metric. This pipeline can be broken down into 3 stages:

1) Compression
* Input: Set of sequences (fasta/fastq files)
* Output: Distance matrix of NCD values

2) Clustering
* Input: Distance matrix of NCD values
* Output: Table of cluster assignments

3) Visualization
* Input: Table of cluster assignments & distance matrix of NCD values
* Output: Sequence similarity network

I currently have a naive implementation of the compression stage, however it requires many optimizations. During development, we will benchmark each stage of the pipeline and compare the runtime and memory requirements to other similarity metrics. We will also apply our pipeline towards datasets of pathogenic bacteria & viruses. Finally, I would like to manage this pipeline into a conda package to make it easy for users to run on their own machines.

## Installation and Dependencies
Dependencies are not yet finalized, and will be updated here. Most likely we will require:
- [`python3`](https://python.org)
- [`numpy`](https://numpy.org)
- [`scipy`](https://scipy.org)
- [`pandas`](https://pandas.pydata.org)
- [`matplotlib`](https://matplotlib.org)
- [`sklearn`](http://scikit-learn.org/stable/)
- [`py-lz4framed`](https://github.com/Iotic-Labs/py-lz4framed)

### Set up a virtual env (optional)
To install virtualenv use the following command in your terminal:

    pip install virtualenv

Then in the directory you want to use, create a virtualenv named env:

    virtualenv -p python3.6 env

And then activate the environment with:

    source env/bin/activate

You can leave the virtualenv at any time with the command:

    deactivate

### Install the dependencies

To install the dependencies:

    pip install -e .
