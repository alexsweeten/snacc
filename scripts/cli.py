import click
import os
import concurrent.futures
from itertools import product
from tqdm import tqdm
import pandas as pd

from pairwise_ncd import return_byte, compressed_size, compute_distance

def compute_distance(comparison, algorithm= "lzma"):
    #Compute a distance between a and b
    sequences = return_byte(open(comparison[0].read()), open(comparison[1].read()))
    sizes = compressed_size(sequences, algorithm=algorithm)
    ncd = compute_distance(sizes[0], sizes[1], sizes[2])
    return (comparison[0], comparison[1], ncd)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-f", "--fasta", type=click.Path(dir_okay=False, exists=True, resolve_path=True), multiple=True, help="FASTA file containing sequence to compare")
@click.option("-d", "--directory", "directories", type=click.Path(dir_okay=True, file_okay=False, exists=True, resolve_path=True), multiple=True, help="Directory containing FASTA files to compare")

@click.option("-N", "--numCPU", "numCPU",default= 1, help="Number of CPU cores available")
@click.option("-o", "--output", type=click.Path(dir_okay=False, exists=False), help="The location for the output CSV file")
def cli(fasta, directories,numCPU, output):


    # generate a list of absolute paths containing the files to be compared
    files = list(fasta)

    for directory in directories:
        for dirpath, _, filenames in os.walk(directory):
            for f in filenames:
                files.append(os.path.abspath(os.path.join(dirpath, f)))
    files = list(set(files)) # remove duplicates
    comparisons = tqdm(list(product(files, repeat=2)))

    executor = concurrent.futures.ProcessPoolExecutor(max_workers=numCPU)
    distances = [res for res in executor.map(compute_distance,comparisons)]

    df = pd.DataFrame(distances, columns=["file", "file2", "ncd"])#.to_csv("out.csv", index=False)

    df.pivot(index='file', columns='file2', values='ncd').to_csv(output)


if __name__ == "__main__":
    cli()
