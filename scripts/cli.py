import click
import os
import concurrent.futures
from itertools import product
from tqdm import tqdm
import pandas as pd
from Bio import SeqIO
from pairwise_ncd import return_byte, compressed_size, compute_distance


def extract_sequences(filepath, reverse_complement=False):
    seq = ""
    for seq_record in SeqIO.parse(filepath, "fasta"):
        if reverse_complement:
            seq += str(seq_record.seq.reverse_complement())
        else:
            seq += str(seq_record.seq)
    return seq


def tqdm_parallel_map(showProgress, executor, fn, *iterables, **kwargs):
    """
    Equivalent to executor.map(fn, *iterables),
        but displays a tqdm-based progress bar.
        Does not support timeout or chunksize as executor.submit is used internally
    **kwargs is passed to tqdm.
    """
    futures_list = []
    for iterable in iterables:
        futures_list += [executor.submit(fn, i) for i in iterable]
    if showProgress:
        print("Processes submitted, starting compression distance calculation...")
        for f in tqdm(concurrent.futures.as_completed(futures_list), total=len(futures_list), **kwargs):
            yield f.result()
    else:
        for f in concurrent.futures.as_completed(futures_list):
            yield f.result()

def compute_parallel(comparison, algorithm, saveCompression="", reverse_complement=False, bwtDiskDict=dict()):
    #Compute a distance between a and b
    sequences = return_byte(extract_sequences(comparison[0], reverse_complement=reverse_complement),
                            extract_sequences(comparison[1], reverse_complement=reverse_complement))
    sizes = compressed_size(sequences, algorithm, saveCompression, comparison)
    ncd = compute_distance(sizes[0], sizes[1], sizes[2])
    return comparison[0], comparison[1], ncd


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-f", "--fasta", type=click.Path(dir_okay=False, exists=True, resolve_path=True), multiple=True, help="FASTA file containing sequence to compare")
@click.option("-d", "--directory", "directories", type=click.Path(dir_okay=True, file_okay=False, exists=True, resolve_path=True), multiple=True, help="Directory containing FASTA files to compare")
@click.option("-n", "--num-threads", "numThreads", type=int, default=None, help="Number of Threads to use (default 5 * number of cores)")
@click.option("-o", "--output", required=True, type=click.Path(dir_okay=False, exists=False), help="The location for the output CSV file")
@click.option("-s", "--save-compression", "saveCompression", type=click.Path(dir_okay=True, file_okay=False, resolve_path=True), help="Save compressed sequence files to the specified directory")
@click.option("-c", "--compression", default="lzma", type=click.Choice(['bwt-disk', 'lzma', 'gzip', 'bzip2', 'zlib', 'lz4', 'snappy']), help="The compression algorithm to use. Defaults to lzma.")
@click.option("-p", "--show-progress", "showProgress", default=True, type=bool, help="Whether to show a progress bar for computing compression distances")
@click.option("-r", "--reverse-complement", is_flag=True, default=False, help="Whether to use the reverse complement of the sequence")
@click.option("-bM", "--bwte-mem", "bwteMem", type=int, help="BWT-Disk: internal memory to be used in bwt-disk in MB (def. 256 MB)")
@click.option("-bC", "--bwte-compress", "bwteCompress", type=click.Choice(['gzip', 'lzma', 'rle-range-encode', 'dna5-symbol']), default='lzma', help="BWT-Disk: compression to be used after running bwt-disk")
def cli(fasta, directories, numThreads, compression, showProgress, saveCompression, output, reverse_complement, bwteMem, bwteCompress):
    bwtDict = {}
    if compression == 'bwt-disk': #construct input dictionary to bwt-disk executable
        filters = {
            '': 0,
            'gzip': 1,
            'rle-range-encode': 2,
            'dna5-symbol': 3,
            'lzma': 4
        }
        bwtDict = {
            'bwt_mem': ['-m', bwteMem],
            'bwt_compress': ['-b', filters[bwteCompress]]
        }
    # generate a list of absolute paths containing the files to be compared
    files = list(fasta)

    for directory in directories:
        for dirpath, _, filenames in os.walk(directory):
            for f in filenames:
                files.append(os.path.abspath(os.path.join(dirpath, f)))
    files = list(set(files)) # remove duplicates
    files = [_file for _file in files if _file.endswith(("fa", "fasta", "faa", "fna"))]
    comparisons = tqdm(list(product(files, repeat=2)))

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=numThreads)
    distances = tqdm_parallel_map(showProgress, executor, lambda x: compute_parallel(x, algorithm=compression,
                                                                                     saveCompression=saveCompression, reverse_complement=reverse_complement, bwtDiskDict=bwtDict), comparisons)

    df = pd.DataFrame(distances, columns=["file", "file2", "ncd"])#.to_csv("out.csv", index=False)

    df.pivot(index='file', columns='file2', values='ncd').to_csv(output)


if __name__ == "__main__":
    cli()
