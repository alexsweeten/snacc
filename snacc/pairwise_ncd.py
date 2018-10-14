import sys
import lzma
import lz4framed
import bz2
import zlib
import os
import subprocess
import bz2
import gzip
import lz4framed
from Bio import SeqIO
from pathlib import Path

def runBwtDisk(filename, inputs, save_directory=None):
	"""
	read in fasta and bwte options to run bwt-disk with compression and filters

		Args:
			filename (pathlib.Path): A path object for the sequence's original fasta
			inputs (dict): a dictionary of options to bwte executable

		Retuns:
			(file): The return statement. A binary file object of the BWT-Disk
					transformed sequence thats compressed
	"""
	if type(filename) == tuple:
		out_file = filename[0].stem + filename[1].name
		seqs = extract_sequences(filename[0]) + extract_sequences(filename[1])
		f = open(out_file,'w+')
		f.write(seqs)
		f.close()
	else:
		out_file = filename.absolute()
	extensions = ['','.gz', '.rrc', '.atn', '.lzma']
	basedir = os.path.abspath(os.path.dirname(__file__))
	bwte_exec = "." + os.path.join(basedir,'bwte')
	cmd = [bwte_exec]
	for key in inputs:
		cmd += inputs[key]
	subprocess.run(cmd + [out_file])

	output = str(out_file) + ".bwt" + extensions[int(inputs['bwte-compress'][1])]
	result = open(output, 'rb').read()
	if type(filename) == tuple:
		if save_directory:
			subprocess.run(['mv','-t',save_directory,output, output + ".aux"])
			subprocess.run(['rm',out_file])
		else:
			subprocess.run(['rm',output, out_file, output + ".aux"])
	else:
		if save_directory:
			subprocess.run(['mv','-t',save_directory,output, output +".aux"])
		else:
			subprocess.run(['rm',output, output + ".aux"])
	return result

def bwt(s):
    """
    Code for Burrows-Wheeler Transformation. The code is adapted from Rosetta Code.
    """
    # uncomment below if you plan on decoding
    #s = "\002" + s + "\003"
    table = sorted(s[i:] + s[:i] for i in range(len(s)))
    last_column = [row[-1:] for row in table]
    return "".join(last_column)


def extract_sequences(filepath, reverse_complement=False, BWT=False):
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
    if BWT:
        return bwt(seq)
    else:
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

    # check if already compressed @TODO
    sequence = bytes(extract_sequences(filename, reverse_complement=reverse_complement, BWT=BWT), encoding="utf-8")
    extension = {
        "bwt-disk": "",
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
    elif algorithm == 'bwt-disk':
        compressed_seq = runBwtDisk(filename, bwte_inputs)

    if save_directory:
        if type(filename) == tuple:
            out_file = filename[0].stem + filename[1].name
        else:
            out_file = filename.name
        with open(os.path.join(save_directory.absolute(), out_file + extension[algorithm]), 'wb') as f:
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
