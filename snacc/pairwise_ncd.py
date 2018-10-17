
import bz2
import gzip
import lzma
import os
import sys
import tempfile
import zlib
import subprocess
import lz4framed
from Bio import SeqIO
from pathlib import Path

def runBwtDisk(seq, inputs, extension):
    """
    read in fasta and bwte options to run bwt-disk with compression and filters

        Args:
            seq (str): A sequence stripped of all headers etc.
            inputs (dict): a dictionary of options to bwte executable

        Retuns:
            (file): The return statement. A binary file object of the BWT-Disk
                    transformed sequence thats compressed

    """
    basedir = os.path.abspath(os.path.dirname(__file__))
    bwte_exec =  os.path.join(basedir,'../bin/bwt_disk/bwte')
    cmd = [bwte_exec]
    for key in inputs:
        cmd += inputs[key]
    with tempfile.NamedTemporaryFile(mode='w+') as f:
        f.write(seq)
        subprocess.run(cmd + [f.name])
        result = open(f.name + extension, "rb").read()
    toRemove = [f.name + extension, f.name + extension + ".aux"]
    subprocess.run(["rm"] + toRemove)
    return result


def extract_sequences(filepath, reverse_complement=False):
    if type(filepath) == tuple:
        return extract_sequences(filepath[0]) + extract_sequences(filepath[1])
    seq = ""
    for seq_record in SeqIO.parse(filepath.absolute(), "fasta"):
        if reverse_complement:
            seq += str(seq_record.seq.reverse_complement())
        else:
            seq += str(seq_record.seq)
    if not seq:
        raise ValueError(f"No sequence extracted. Ensure that file {filepath.absolute()} contains a proper FASTA definition line (i.e. a line that starts with '>sequence_name').")
    return seq


def compressed_size(filename, algorithm, reverse_complement=False, save_directory=None, BWT=False, bwte_inputs = {}):
    '''

    Args:
        filename (pathlib.Path)
        algorithm (str)
                reverse_complement(bool, optional)
        save_directory (pathlib.Path, optional)

    Returns
        (pathlib.Path,int): the number of bytes in the compressed file
    '''

    sequence = extract_sequences(filename, reverse_complement=reverse_complement)
    extension = {
        "lzma": ".lzma",
        "gzip": ".gz",
        "bzip2": ".bz2",
        "zlib": ".ZLIB",
        "lz4": ".lz4",
        "bwt-disk-rle-range": ".bwt.rrc",
        "bwt-disk-dna5-symbol": ".bwt.atn"
    }
    file_ext = extension[algorithm]

    if BWT:
        if "bwt" not in algorithm:
            file_ext = ".bwt" + file_ext
            sequence = runBwtDisk(sequence,bwte_inputs, ".bwt")
    else:
        if "bwt" not in algorithm:
            sequence = bytes(sequence, encoding = "utf-8")

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
    elif algorithm == 'bwt-disk-rle-range':
        compressed_seq = runBwtDisk(sequence, bwte_inputs, file_ext)
    elif algorithm == 'bwt-disk-dna5-symbol':
        compressed_seq = runBwtDisk(sequence, bwte_inputs, file_ext)

    if save_directory:
        if type(filename) == tuple:
            out_file = filename[0].stem + filename[1].name
        else:
            out_file = filename.name
        with open(os.path.join(save_directory.absolute(), out_file + file_ext), 'wb') as f:
            f.write(compressed_seq)

    return (filename, sys.getsizeof(compressed_seq))


#calculates NCD for 2 sequence sizes and their concatenation size
def compute_distance(x, y, cxy, cyx):
    if x > y:
        return min((cxy - y) / x, (cyx - y) / x)
    elif y > x:
        return min((cxy - x) / y, (cyx - x) / y)
    else:
        return min((cxy - x) / x, (cyx - x) / x)
