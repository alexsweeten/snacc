import sys
import lzma
import lz4framed
import bz2
import zlib
import getopt
from sklearn.cluster import AgglomerativeClustering
import sys
import argparse
import os
import bz2
import gzip
import snappy
import lz4framed

def bwt(s):
    """
    Code for Burrows-Wheeler Transformation. The code is adapted from Rosetta Code.
    """
    s = "\002" + s + "\003"
    table = sorted(s[i:] + s[:i] for i in range(len(s)))
    last_column = [row[-1:] for row in table]
    return "".join(last_column)


#given 2 sequence strings, returns sequences + concatenation as an object of bytes
def return_byte(sequence1, sequence2):
    seq1 = bytes(sequence1, 'utf-8')
    seq2 = bytes(sequence2, 'utf-8')
    return seq1, seq2, seq1 + seq2

def extract_sequences(filepath, reverse_complement=False):
	if type(filepath) == tuple:
		return extract_sequences(filepath[0]) + extract_sequences(filepath[1])
    seq = ""
    for seq_record in SeqIO.parse(filepath.absolute(), "fasta"):
        if reverse_complement:
            seq += str(seq_record.seq.reverse_complement())
        else:
            seq += str(seq_record.seq)
    return seq


def compressed_size(filename, algorithm, reverse_complement=False, save_directory=None):
    '''

    Args:
        filename (pathlib.Path)
        algorithm (str)
		reverse_complement(bool, optional)
        save_directory (pathlib.Path, optional)

    Returns
        (pathlib.Path,int): the number of bytes in the compressed file
    '''

    # check if already compressed @TODO
	sequence = extract_sequences(filename, reverse_complement=reverse_complement)
    extension = {
        "lzma": ".lzma",
        "gzip": ".gz",
        "bzip2": ".bz2",
        "zlib": ".ZLIB",
        "lz4": ".lz4"
    }
    if algorithm == "lzma":
        compressed_seq = lzma.compress(sequence)
    elif algorithm == "gzip":
        compressed_seq = gzip.compress(sequence)
    elif algorithm == "bzip2":
        compressed_seq = bz2.compress(sequence)
    elif algorithm == "zlib":
        compressed_seq = zlib.compress(sequence)
    elif algorithm == "lz4":
        compressed_seq = lz4framed.compress(sequence)
    elif algorithm == 'snappy':
        compressed_seq = snappy.compress(sequence)

    if save_directory:
        with open(os.path.join(save_directory.absolute(), filename.name + extension[algorithm]), 'wb') as f:
            f.write(compressed_seq)

    return (filename,sys.getsizeof(compressed_seq))
#calculates NCD for 2 sequence sizes and their concatenation size
def compute_distance(x, y, cxy, cyx):
    if x > y:
        return min((cxy - y) / x, (cyx - y) / x)
    elif y > x:
        return min((cxy - x) / y, (cyx - x) / y)
    else:
        return min((cxy - x) / x, (cyx - x) / x)
