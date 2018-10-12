import click
from Bio.SeqIO import parse
from pprint import pprint

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument("FASTA", type=click.Path(dir_okay=False, exists=True, resolve_path=True), nargs=-1)
def cli(fasta):
    # for now, just print the name of the files we got
    pprint(fasta)


if __name__ == "__main__":
    cli()
