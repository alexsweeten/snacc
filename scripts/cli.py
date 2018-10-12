import click

from pairwise_ncd import return_byte, compressed_size, compute_distance

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument("FASTA", type=click.Path(dir_okay=False, exists=True, resolve_path=True), nargs=2)
def cli(fasta):

    #convert input sequences into bytes
    sequences = return_byte(open(fasta[0]).read(), open(fasta[1]).read())

    #compress input sequences
    sizes = compressed_size(sequences)

    #compute ncd values
    ncd = compute_distance(sizes[0], sizes[1], sizes[2])

    print(ncd)


if __name__ == "__main__":
    cli()
