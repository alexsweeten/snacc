"""
computeGenCompress.py:
    Function for compressing a raw fasta file using the gencompress algorithm
"""
import os
import subprocess


def GenCompress(fasta):
    """
    Run GenCompress executable on a fasta file and return a fileobjet containing
    the compressed format

    Args:
        fasta(str): A file path to the fasta file to compress

    Returns:
        (file): The return statement. A compressed fasta file object
    """
    (fname, ext) = os.path.splitext(fasta)
    out_name = ".GEN"
    executable =  "../bin/GenCompress"
    cmd = "./" + executable
    subprocess.run([cmd, fasta])
    output = open(out_name)
    subprocess.run(["rm", out_name])
    return output
