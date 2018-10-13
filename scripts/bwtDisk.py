"""
bwtDisk.py:
	backend to bwt-disk binary available at https://people.unipmn.it/manzini/bwtdisk/ .
	This manages the bwt-disk options separate from the other compression algorithms
	and protocols
"""
import os
import subprocess

def runBwtDisk(fasta, **kwargs):
	cmd = ["./../bin/bwt_disk/bwte"]
	for key in kwargs:
		cmd += kwargs[key]
	subprocess.run(cmd)
	output = kwargs['bwt_out'] + extensions[kwargs['bwt_compress']]
	return open(output, 'rb')
