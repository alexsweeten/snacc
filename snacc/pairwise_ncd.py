
import bz2
import gzip
import lzma
import os
import subprocess
import sys
import tempfile
import zlib
from pathlib import Path

import lz4framed
from Bio import SeqIO

def extract_sequences(sequences, reverse_complement=False):
    """
    Extracts and concatenates the sequences within FASTA files.

    Args:
        sequences (pathlib.Path or tuple): Either the :obj:`~pathlib.Path` of the FASTA file or a tuple of FASTA files from which to extract sequences.
        reverse_complement (bool, optional): Whether to take the reverse complement of the sequences.

    Returns:
        str: The concatenated sequences within the FASTA file(s).

    Raises:
        ValueError: When the FASTA file is malformed.
    """
    if type(sequences) == tuple:
        return extract_sequences(sequences[0], reverse_complement=reverse_complement) + extract_sequences(sequences[1], reverse_complement=reverse_complement)
    seq = ""
    for seq_record in SeqIO.parse(sequences.absolute(), "fasta"):
        if reverse_complement:
            seq += str(seq_record.seq.reverse_complement())
        else:
            seq += str(seq_record.seq)
    if not seq:
        raise ValueError(f"No sequence extracted. Ensure that file {sequences.absolute()} contains a proper FASTA definition line (i.e. a line that starts with '>sequence_name').")
    return seq


def compressed_size(sequences, algorithm, reverse_complement=False, save_directory=None, BWT=False, bwte_inputs={}):
    '''
    Calculates the compressed size of the sequences in a file or tuple of files.

    Args:
        sequences (pathlib.Path or tuple): Either the :obj:`~pathlib.Path` of the FASTA file to compress or a tuple of FASTA files to concatenate and compress.
        algorithm (str): Which algorithm to compress the file with. Valid options are [``lzma``, ``gzip``, ``bzip2``, ``zlib``, ``lz4``, ``bwt-disk-rle-range``, ``bwt-disk-dna5-symbol``].
        reverse_complement(bool, optional): Whether to take the reverse complement of the sequences in the file.
        save_directory (pathlib.Path, optional): If given, where to save the compressed file.

    Note:
        The entire file is not compressed, just the sequences within it.

    Returns:
        tuple: A tuple whose zeroth element is ``sequences`` (_i.e._ either a :obj:`~pathlib.Path` or a tuple) and first element is the number of bytes in the compressed file.
    '''

    sequence = extract_sequences(sequences, reverse_complement=reverse_complement)
    extension = {
        "lzma": ".lzma",
        "gzip": ".gz",
        "bzip2": ".bz2",
        "zlib": ".ZLIB",
        "lz4": ".lz4"
    }
    file_ext = extension[algorithm]

    sequence = bytes(sequence, encoding="utf-8")

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

    if save_directory:
        if type(sequences) == tuple:
            out_file = sequences[0].stem + sequences[1].name
        else:
            out_file = sequences.name
        with open(os.path.join(save_directory.absolute(), out_file + file_ext), 'wb') as f:
            f.write(compressed_seq)

    return (sequences, sys.getsizeof(compressed_seq))


def compute_distance(x, y, cxy, cyx):
    """
    Calculates normalized compression distance for two files and their concatenations.

    Args:
        x (int): The size (in bytes) of the first file.
        y (int): The size (in bytes) of the second file.
        cxy (int): The size (in bytes) of the second file concatenated onto the first.
        cyx (int): The size (in bytes) of the first file concatenated onto the second.

    Returns:
        int: The distance between the two files.
    """
    if x > y:
        return min((cxy - y) / x, (cyx - y) / x)
    elif y > x:
        return min((cxy - x) / y, (cyx - x) / y)
    else:
        return min((cxy - x) / x, (cyx - x) / x)
