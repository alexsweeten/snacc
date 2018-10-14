import concurrent.futures
import itertools
import os
import sys
from datetime import datetime
from pathlib import Path
import webbrowser

import click
import jinja2
import lz4framed
import pandas as pd
import sklearn
import umap
from markdown import markdown
from tqdm import tqdm

from .pairwise_ncd import compressed_size, compute_distance
from .version import __version__


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-f", "--fasta", type=click.Path(dir_okay=False, exists=True, resolve_path=True), multiple=True, help="FASTA file containing sequence to compare.")
@click.option("-d", "--directory", "directories", type=click.Path(dir_okay=True, file_okay=False, exists=True, resolve_path=True), multiple=True, help="Directory containing FASTA files to compare.")
@click.option("-n", "--num-threads", "numThreads", type=int, default=None, help="Number of Threads to use (default 5 * number of cores).")
@click.option("-o", "--output", type=click.Path(dir_okay=False, exists=False), help="The location for the output CSV file.")
@click.option("-s", "--save-compression", "saveCompression", type=click.Path(dir_okay=True, file_okay=False, resolve_path=True), default=None, help="Save compressed sequence files to the specified directory.")
@click.option("-c", "--compression", default="lzma", type=click.Choice(['lzma', 'gzip', 'bzip2', 'zlib', 'lz4', 'snappy', 'bwt-disk']), help="The compression algorithm to use. Defaults to lzma.")
@click.option("--show-progress/--no-show-progress", "showProgress", default=True, help="Whether to show a progress bar for computing compression distances.")
@click.option("-r", "--reverse_complement", is_flag=True, default=False, help="Whether to use the reverse complement of the sequence.")
@click.option("-b", "--burrows-wheeler", "BWT", is_flag=True, default=False, help="Whether to compute the Burrows-Wheeler Tranform prior to compression and reverse complement (default 256 MB).")
@click.option("-bM", "--bwte-mem", "bwteMem", type=int, default=256, help="BWT-Disk option: The amount of memory in MB for use in the bwt-disk executable.")
@click.option("-bC", "--bwte-compress", "bwteCompress", type=click.Choice(['None', 'gzip', 'rle-range-encoding', 'dna5-symbol', 'lzma']), default='gzip', help="BWT-Disk Option: The compression to use when calling bwt-disk before compression, may require separate libraries if not using default.")
@click.option("-l", "--log-type", default="html", type=click.Choice(["html", "md"]), help="The output format for the report.")
@click.option("--no-show", default=False, is_flag=True, help="If mode is html, use this flag to prevent automatically opening the log in the browser.")
def cli(fasta, directories, numThreads, compression, showProgress, saveCompression, output, reverse_complement, BWT, bwteMem, bwteCompress, log_type, no_show):
    start_time = datetime.now()

    if saveCompression:
        saveCompression = Path(saveCompression)
    # Map bwte inputs to options for the bwte executabl
    compressions = {
        'None': 0,
        'gzip': 1,
        'rle-range-encoding': 2,
        'dna5-symbol': 3,
        'lzma': 4
    }
    bwte_inputs = {
        'bwte-mem': ['-m', str(bwteMem)],
        'bwte-compress': ['-b', str(compressions[bwteCompress])]
    }

    output = Path(output) # make the output into a Path object

    # generate a list of absolute paths containing the files to be compared
    files = [Path(f) for f in fasta]

    # get all the files in the passed directories
    for directory in directories:
        directory = Path(directory)
        for f in directory.iterdir():
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
                                             BWT=BWT,
                                             bwte_inputs=bwte_inputs),
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
                                                   BWT=BWT,
                                                   bwte_inputs=bwte_inputs),
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
                                           table=df.to_html(index=False),
                                           duration=datetime.now() - start_time,
                                           output_path=output.absolute())

    if log_type == "html":
        rendered = markdown(rendered)
        rendered = jinja2.Template(rendered_html).render(body=rendered)

    print(rendered, file=open(output.stem + "." + log_type, "w"))
    if not no_show and log_type == "html":
        webbrowser.open("file://" + str(output.parent.absolute()) + "/" + output.stem + "." + log_type)


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


rendered_html = '''
<!doctype html>
<html lang="en">
  <head>
    <!-- Required meta tags -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">

    <!-- Bootstrap CSS -->
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css" integrity="sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO" crossorigin="anonymous">
    <title>snacc analysis</title>
  </head>
  <body>
  {{ body }}
  </body>
<style>
body {
    padding-left: 50px;
    padding-top: 25px;
}
</style>
</html>
'''
log = '''#snacc analysis
* Analysis time: {{time}}
* Analysis duration: {{duration}}
* Compression method: {{method}}
* Reverse complement: {{rev_comp}}
* Burrows-Wheeler transform: {{bwt}}
* Output filepath: {{output_path}}

### Analyzed Files
{% for _file in files -%}
* {{_file}}
{% endfor %}

### Distance Matrix
{{table}}

### Version Information
* Python: {{py_version}}
* snacc: {{snacc_version}}
* scikit-learn: {{sklearn_version}}
* py-lz4framed: {{lz4framed_version}}
* umap-learn: {{umap_version}}
'''

# Utility function
def convert_html_to_pdf(sourceHtml, outputFilename):
    # open output file for writing (truncated binary)
    resultFile = open(outputFilename, "w+b")

    # convert HTML to PDF
    pisaStatus = pisa.CreatePDF(sourceHtml,      # the HTML to convert
                                dest=resultFile) # file handle to recieve result
    resultFile.close()


if __name__ == "__main__":
    cli()
