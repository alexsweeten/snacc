import concurrent.futures
import itertools
import os
import sys
from datetime import datetime
from pathlib import Path

import click
import jinja2
import lz4framed
import pandas as pd
import sklearn
import umap
from markdown import markdown
from tqdm import tqdm
from xhtml2pdf import pisa

from .pairwise_ncd import compressed_size, compute_distance
from .version import __version__


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-f", "--fasta", type=click.Path(dir_okay=False, exists=True, resolve_path=True), multiple=True, help="FASTA file containing sequence to compare.")
@click.option("-d", "--directory", "directories", type=click.Path(dir_okay=True, file_okay=False, exists=True, resolve_path=True), multiple=True, help="Directory containing FASTA files to compare.")
@click.option("-n", "--num-threads", "numThreads", type=int, default=None, help="Number of Threads to use (defaults to 5 * number of cores).")
@click.option("-o", "--output", type=click.Path(dir_okay=False, exists=False), help="The location for the output CSV file.")
@click.option("-s", "--save-compression", "saveCompression", type=click.Path(dir_okay=True, file_okay=False, resolve_path=True), default=None, help="Save compressed sequence files to the specified directory.")
@click.option("-c", "--compression", default="lzma", type=click.Choice(['lzma', 'gzip', 'bzip2', 'zlib', 'lz4']), help="The compression algorithm to use. Defaults to lzma.")
@click.option("-p", "--show-progress", "showProgress", default=True, type=bool, help="Whether to show a progress bar for computing compression distances.")
@click.option("-r", "--reverse_complement", is_flag=True, default=False, help="Whether to use the reverse complement of the sequence.")
@click.option("-b", "--burrows-wheeler", "BWT", is_flag=True, default=False, help="Whether to compute the Burrows-Wheeler Transform prior to compression and reverse complement.")
@click.option("-l", "--log-mode", type=click.Choice(["html", "md"]), help="The output format for the report.")
def cli(fasta, directories, numThreads, compression, showProgress, saveCompression, output, reverse_complement, BWT, log_mode):
    if saveCompression:
        saveCompression = Path(saveCompression)
    # generate a list of absolute paths containing the files to be compared
    files = [Path(f) for f in fasta]

    for directory in directories:
        for dirpath, _, filenames in os.walk(directory):
            for f in filenames:
                f = Path(os.path.join(dirpath, f))
                if f.suffix.lower() in [".fasta", ".fna", ".fa", ".faa"]:
                    files.append(f)
    files = list(set(files)) # remove any duplicates

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=numThreads)

    #compute compressed sizes of individual sequences
    print("Compressing individual files...")
    compressed_sizes = tqdm_parallel_map(executor,
                                         lambda x: compressed_size(
                                             filename=x,
                                             algorithm=compression,
                                             save_directory=saveCompression,
                                             reverse_complement=reverse_complement,
                                             BWT=BWT),
                                         showProgress,
                                         files)
    compressed_dict = dict(compressed_sizes) # {PATH: compressed size}

    # compute compressed sizes of all ordered pairs of sequences
    print("Compressing pairs...")
    compressed_pairs_sizes = tqdm_parallel_map(executor,
                                               lambda x: compressed_size(
                                                   filename=x,
                                                   algorithm=compression,
                                                   save_directory=saveCompression,
                                                   reverse_complement=reverse_complement,
                                                   BWT=BWT),
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

    rendered = jinja2.Template(log).render(time=datetime.now(),
                                           method=compression,
                                           py_version=str(sys.version.replace("\n", "")),
                                           snacc_version=__version__,
                                           umap_version=umap.__version__,
                                           sklearn_version=sklearn.__version__,
                                           lz4framed_version=lz4framed.__version__,
                                           files=[str(_file.absolute()) for _file in files],
                                           bwt=BWT,
                                           rev_comp=reverse_complement,
                                           table=df.to_html(index=False))

    print(rendered, file=open("test.md", "w"))

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


log = '''# `snacc` Analysis
* Analysis time: {{time}}
* Compression method: {{method}}
* Reverse complement: {{rev_comp}}
* Burrows-Wheeler transform: {{bwt}}

## Analyzed Files
{% for _file in files -%}
* {{_file}}
{% endfor %}

## Result
{{table}}

## Version Information
* Python: {{py_version}}
* snacc: {{snacc_version}}
* scikit-learn: {{sklearn_version}}
* py-lz4framed: {{lz4framed_version}}
* umap-learn: {{umap_version}}

'''


if __name__ == "__main__":
    cli()
