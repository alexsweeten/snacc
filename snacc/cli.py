import concurrent.futures
import itertools
import os
import sys
from datetime import datetime
from pathlib import Path

import click
import lz4framed
import pandas as pd
from tqdm import tqdm

from .pairwise_ncd import compressed_size, compute_distance
from .version import __version__


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument("sequences",
                type=click.Path(exists=True,
                                resolve_path=True),
                nargs=-1)
@click.option("-f", "--fasta",
              type=click.Path(dir_okay=False,
                              exists=True,
                              resolve_path=True),
              multiple=True,
              hidden=True,
              help="FASTA file containing sequence to compare.")
@click.option("-d", "--directory", "directories",
              type=click.Path(dir_okay=True,
                              file_okay=False,
                              exists=True,
                              resolve_path=True),
              multiple=True,
              hidden=True,
              help="Directory containing FASTA files to compare.")
@click.option("-n", "--num-threads", "numThreads",
              type=int,
              default=None,
              help="Number of Threads to use (default 5 * number of cores).")
@click.option("-o", "--output",
              type=click.Path(dir_okay=False,
                              exists=False),
              help="The location for the output CSV file.",
              prompt="Output CSV path")
@click.option("-s", "--save-compression", "saveCompression",
              type=click.Path(dir_okay=True,
                              file_okay=False,
                              resolve_path=True),
              default=None,
              help="Save compressed sequence files to the specified directory.")
@click.option("-c", "--compression",
              default="lzma",
              type=click.Choice(['lzma',
                                 'gzip',
                                 'bzip2',
                                 'zlib',
                                 'lz4',
                                 'bwt-disk-rle-range',
                                 'bwt-disk-dna5-symbol']),
              help="The compression algorithm to use. Defaults to lzma.")
@click.option("--show-progress/--no-show-progress", "showProgress",
              default=True,
              help="Whether to show a progress bar for computing compression distances.")
@click.option("-r", "--reverse_complement",
              is_flag=True,
              default=False,
              help="Whether to use the reverse complement of the sequence.")
@click.option("--log/--no-log", "log",
              default=True,
              help="Whether to save a log.")
def cli(sequences,
        fasta,
        directories,
        numThreads,
        compression,
        showProgress,
        saveCompression,
        output,
        reverse_complement,
        log_type,
        show_log):
    start_time = datetime.now()

    if fasta or directories:
        click.secho("Warning: the -f and -d flags are deprecated and will be removed before release. Please pass files and paths directly without the flags.", fg="yellow")

    if saveCompression:
        saveCompression = Path(saveCompression)
    
    output = Path(output) # make the output into a Path object

    # generate a list of absolute paths containing the files to be compared
    sequences = [Path(sequence) for sequence in sequences]
    files = [path for path in sequences if path.is_file()]

    # temporary fix for #80 (to be removed when the -d and -f flags are removed)
    files.extend([Path(_f) for _f in fasta])
    sequences.extend([Path(_f) for _f in directories])

    # get all the files in the passed directories
    for directory in [path for path in sequences if path.is_dir()]:
        for f in directory.iterdir():
            if f.suffix.lower() in [".fasta", ".fna", ".fa", ".faa", ".fsa"]:
                files.append(f)
    files = sorted(list(set(files)), key=lambda x: str(x.absolute())) # remove any duplicates and sort for cleaner log output

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=numThreads)

    #compute compressed sizes of individual sequences
    click.secho("Compressing individual files...", fg="green")
    compressed_sizes = tqdm_parallel_map(executor,
                                         lambda x: compressed_size(
                                             sequences=x,
                                             algorithm=compression,
                                             save_directory=saveCompression,
                                             reverse_complement=reverse_complement),
                                         showProgress,
                                         files)
    compressed_dict = dict(compressed_sizes) # {PATH: compressed size}

    # compute compressed sizes of all ordered pairs of sequences
    click.secho("Compressing pairs...", fg="green")
    compressed_pairs_sizes = tqdm_parallel_map(executor,
                                               lambda x: compressed_size(
                                                   sequences=x,
                                                   algorithm=compression,
                                                   save_directory=saveCompression,
                                                   reverse_complement=reverse_complement),
                                               showProgress,
                                               itertools.product(compressed_dict.keys(), repeat=2))

    compressed_pairs_dict = dict(compressed_pairs_sizes) # {(A, B): size, (B, A): size,...}

    distances = {}
    for pair in itertools.product(compressed_dict.keys(), repeat=2):
        distances[pair] = compute_distance(compressed_dict[pair[0]],
                                           compressed_dict[pair[1]],
                                           compressed_pairs_dict[(pair[0], pair[1])],
                                           compressed_pairs_dict[(pair[1], pair[0])])

    distances = list(distances.items())
    distances = [(distance[0][0], distance[0][1], distance[1]) for distance in distances]
    df = pd.DataFrame(distances, columns=["file", "file2", "ncd"])#.to_csv("out.csv", index=False)
    df = df.pivot(index='file', columns='file2', values='ncd')
    df.to_csv(output)
    df = pd.read_csv(output)
    df.columns = list(map(lambda x: Path(x).name, df.columns))
    df.file = df.file.apply(lambda x: Path(x).name)

    if log:
        rendered = log_template.format(time=datetime.now(),
                                       method=compression,
                                       py_version=str(sys.version.replace("\n", "")),
                                       snacc_version=__version__,
                                       umap_version=umap.__version__,
                                       sklearn_version=sklearn.__version__,
                                       lz4framed_version=lz4framed.__version__,
                                       bwt=BWT,
                                       rev_comp=reverse_complement,
                                       duration=datetime.now() - start_time,
                                       output_path=output.absolute())

        with open(output.stem + ".md", "w") as f:
            print(rendered, file=f)
            for _f in [str(_file.absolute()) for _file in files]:
                print("*", _f, file=f)

def tqdm_parallel_map(executor, fn, showProgress, *iterables, **kwargs):
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
        for f in tqdm(concurrent.futures.as_completed(futures_list), total=len(futures_list), **kwargs):
            yield f.result()
    else:
        for f in concurrent.futures.as_completed(futures_list):
            yield f.result()


log_template = '''# `snacc` Analysis
## Run Information
* Analysis time: {time}
* Analysis duration: {duration}
* Compression method: {method}
* Reverse complement: {rev_comp}
* Output filepath: {output_path}

## Version Information
* Python: {py_version}
* snacc: {snacc_version}
* scikit-learn: {sklearn_version}
* py-lz4framed: {lz4framed_version}
* umap-learn: {umap_version}

## Analyzed Files
'''


if __name__ == "__main__":
    cli()
