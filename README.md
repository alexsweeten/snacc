![logo](https://github.com/SweetiePi/bioncd-hackseq/blob/master/logo/snacc-header.jpg)
# snacc: alignment-free genome comparison utilizing the normalized compression distance
`snacc` (sequence non-alignment compression & comparison) is a program implementing the normalized compression distance (NCD) specifically for biological data. These distances can be used for clustering, or to rapidly infer phylogenies for large sets of genomes.

## Dependencies

- [`python3`](https://python.org)
- [`numpy`](https://numpy.org)
- [`scipy`](https://scipy.org)
- [`pandas`](https://pandas.pydata.org)
- [`matplotlib`](https://matplotlib.org)
- [`sklearn`](http://scikit-learn.org/stable/)
- [`py-lz4framed`](https://github.com/Iotic-Labs/py-lz4framed)
- [`click`](https://click.palletsprojects.com/en/7.x/)
- [`tqdm`](https://pypi.org/project/tqdm/)
- [`biopython`](https://biopython.org/)

## Installation

To install `snacc` directly you may use:
    
```bash
virtualenv env --python=python3.6 # optional, but recommended to create a clean environment
source env/bin/activate # if using virtualenv, activate the env
pip install git+https://github.com/SweetiePi/snacc
```
    
We recommend you create a conda environment for `snacc` and install PathOGiST through conda.
`snacc` requires Python 3.6, so create a conda environment with the right Python version:
```bash
conda create --name snacc python=3.6
```
And then activate the environment and install `snacc`:
```bash
source activate snacc
conda install -c asweeten snacc
```
When inside the `snacc` conda environment, you can verify correct isntallation by running `snacc -h`.

## Examples

1) Most basic usage
```bash
snacc [folder with sequences] -o [output name]
```
2) Intermediate: customize number of threads and compression algorithm
```
snacc -d [folder with sequences] -o [output name] -n 24 -c gzip
```
3) Full control
```bash
snacc \
--directory [folder with sequences] \
--output [output name] \
--num-threads 24 \
--compression lz4 \
--fast-mode True \
--reverse-compliment False
```

## Output Example
### snacc analysis
* Analysis time: 2018-10-14 15:18:17.257619
* Analysis duration: 0:00:26.383997
* Compression method: lz4
* Reverse complement: False
* Burrows-Wheeler transform: False
* Output filepath: test.csv

#### Version Information
* Python: 3.6.0 (v3.6.0:41df79263a11, Dec 22 2016, 17:23:13) [GCC 4.2.1 (Apple Inc. build 5666) (dot 3)]
* snacc: 0.0.1
* scikit-learn: 0.20.0
* py-lz4framed: 0.12.0
* umap-learn: 0.3.5

#### Analyzed Files
* /test_dataset/mysteryGenome_1.fasta
* /test_dataset/mysteryGenome_2.fasta
