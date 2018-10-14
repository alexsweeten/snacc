![logo](https://github.com/SweetiePi/bioncd-hackseq/blob/master/logo/snacc-header.jpg)
# snacc: Compress and compare pathogen genomes without sequence alignment
snacc is a pipeline that implements the normalized compression distance(NCD) specifically for biological data. The workflow primarily consists of 3 stages: compression, clustering, and visualization. The goal of this project is to provide a faster method of comparing large-scale pathogen genomes to conventional alignment-based methods such as BLAST by exploiting the inherent redundancies of the genetic code.

## Workflow

1) Compression
* Input: Set of sequences (fasta/fastq files)
* Output: Distance matrix of NCD values

2) Clustering
* Input: Distance matrix of NCD values
* Output: Table of cluster assignments

3) Visualization
* Input: Table of cluster assignments & distance matrix of NCD values
* Output: Sequence similarity network

## Installation and Dependencies

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
    
## Usage
    $ snacc
    -f --fastatype              :FASTA file containing sequence to compare
    -d --directory              :Directory containing FASTA files to compare
    -n --num-threads            :Number of Threads to use
    -o --output:                :Location for the output CSV file
    -s --save-compression       :(default=None) Save compressed sequence files to the specified directory
    -c --compression            :(default="lzma") The compression algorithm to use. Choose from 'lzma', 'gzip', 'bzip2', 'zlib', 'lz4', and 'snappy'
    -p --show-progress:         :(default=True) Whether to show a progress bar for computing compression distances
    -r --reverse_complement     :(default=False) Whether to use the reverse complement of the sequence
    -b --burrows-wheeler        :(default=False) Whether to compute the Burrows-Wheeler Tranform prior to compression and reverse complement

## Examples
