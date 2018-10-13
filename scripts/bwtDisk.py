"""
bwtDisk.py:
	backend to bwt-disk binary available at https://people.unipmn.it/manzini/bwtdisk/ .
	This manages the bwt-disk options separate from the other compression algorithms
	and protocols
"""
import os
import subprocess

def runBwtDisk(fasta, save, inputs):
	"""
	read in fasta and bwte options to run bwt-disk with compression and filters

		Args:
			fasta (str): The path to a fasta file
			inputs (dict): a dictionary of options to bwte executable

		Retuns:
			(file): The return statement. A binary file object of the BWT-Disk
					transformed sequence thats compressed
	""""
	extensions = ['','.gz', '.rrc', '.atn', '.lzma']
	cmd = ["./../bin/bwt_disk/bwte"]
	for key in inputs:
		cmd += inputs[key]
	subprocess.run(cmd + [fasta])
	output = fasta + ".bwt" + extensions[inputs['bwt_compress']]]
	if
	return open(output, 'rb')
