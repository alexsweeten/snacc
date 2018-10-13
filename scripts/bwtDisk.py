"""
bwtDisk.py:
	backend to bwt-disk binary available at https://people.unipmn.it/manzini/bwtdisk/ .
	This manages the bwt-disk options separate from the other compression algorithms
	and protocols
"""
import os
import subprocess

def runBwtDisk(fasta, **kwargs):
	"""
	read in fasta and bwte options to run bwt-disk with compression and filters

		Args:
			fasta (str): The path to a fasta file
			kwargs (dict): a dictionary of options to bwte executable

		Retuns:
			(file): The return statement. A binary file object of the BWT-Disk
					transformed sequence thats compressed
	""""
	extensions = ['','.gz', '.rrc', '.atn', '.lzma']
	cmd = ["./../bin/bwt_disk/bwte"]
	for key in kwargs:
		cmd += kwargs[key]
	subprocess.run(cmd + [fasta])
	output = fasta + ".bwt" + extensions[kwargs['bwt_compress']]]
	subprocess.run(['mv', kwargs['bwt_out'], output])
	return open(output, 'rb')
