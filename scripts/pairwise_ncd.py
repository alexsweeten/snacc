import sys
import lzma
import bz2
import zlib
import getopt
from sklearn.cluster import AgglomerativeClustering
import sys
import argparse
import os
import bz2
import gzip


def bwt(s: str) -> str:
    """
    Code for Burrows-Wheeler Transformation. The code is adapted from Rosetta Code.
    """
    s = "\002" + s + "\003"
    table = sorted((s[i:] + s[:i] for i in range(len(s))))
    last_column = [row[-1:] for row in table]
    return "".join(last_column)


# given 2 sequence strings, returns sequences + concatenation as an object of bytes
def return_byte(sequence1: str, sequence2: str) -> tuple:
    seq1 = bytes(sequence1, 'utf-8')
    seq2 = bytes(sequence2, 'utf-8')
    return seq1, seq2, seq1 + seq2


def compressed_size(sequences: tuple, algorithm: str, saveCompression: str, comparison: tuple) -> tuple:
	extension = {
		"lzma": ".lzma",
		"gzip": ".gz",
		"bzip2": ".bz2",
		"zlib": ".ZLIB",
		"lz4": ".lz4"
	}
	if algorithm == "lzma":
		compressed_seq1 = lzma.compress(sequences[0])
		compressed_seq2 = lzma.compress(sequences[1])
		compressed_seqconcat = lzma.compress(sequences[2])

	elif algorithm == "gzip":
		compressed_seq1 = gzip.compress(sequences[0])
		compressed_seq2 = gzip.compress(sequences[1])
		compressed_seqconcat = gzip.compress(sequences[2])

	elif algorithm == "bzip2":
		compressed_seq1 = bz2.compress(sequences[0])
		compressed_seq2 = bz2.compress(sequences[1])
		compressed_seqconcat = bz2.compress(sequences[2])

	elif algorithm == "zlib":
		compressed_seq1 = zlib.compress(sequences[0])
		compressed_seq2 = zlib.compress(sequences[1])
		compressed_seqconcat = zlib.compress(sequences[2])

	elif algorithm == "lz4":
		compressed_seq1 = lz4framed.compress(sequences[0])
		compressed_seq2 = lz4framed.compress(sequences[1])
		compressed_seqconcat = lz4framed.compress(sequences[2])

	if saveCompression == "":
		f = open(os.path.join(saveCompression,comparison[0] + extension[algorithm]))
		f.write(compressed_seq1)
		f.close()
		f = open(os.path.join(saveCompression, comparison[0] + extension[algorithm]))
		f.write(compressed_seq2)
		f.close()
		f = open(os.path.join(saveCompression,comparison[2] + extension[algorithm]))
		f.write(compressed_seqconcat)

	compressed_seq1_size = sys.getsizeof(compressed_seq1)
	compressed_seq2_size = sys.getsizeof(compressed_seq2)
	compressed_seqconcat_size = sys.getsizeof(compressed_seqconcat)

	return compressed_seq1_size, compressed_seq2_size, compressed_seqconcat_size

# calculates NCD for 2 sequence sizes and their concatenation size
def compute_distance(x: float, y: float, cxy: float) -> float:
    if x > y:
        distance = ((cxy - y) / x)
    elif y > x:
        distance = ((cxy - x) / y)
    else:
        distance = ((cxy - x) / x)
    return distance
