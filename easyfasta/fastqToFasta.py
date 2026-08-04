"""Command-line entry point converting a FASTQ file into a wrapped FASTA file."""
import argparse
import typing
from collections.abc import Iterable
from typing import TextIO, Generator
from .common import fastq_iter, wrap_sequence


def fastq_to_fasta(fastq: str, fasta: str):
    """
    convert a FASTQ file into a wrapped, multiline FASTA file.

    :param str fastq: the FASTQ file to read
    :param str fasta: path to write the converted fasta to, overwritten if it exists
    :return None: writes the converted fasta to fasta
    """
    with open(fastq) as fi, open(fasta, "w") as fo:
        for record in fastq_iter(fi):
            fo.write(">{}\n{}\n".format(record[0], wrap_sequence(record[1])))

def main():
    """
    parse command-line arguments and run fastq_to_fasta.

    :return None: writes the converted fasta to the path given by --fasta
    """
    parse = argparse.ArgumentParser(description="make a fasta from a fastq")
    parse.add_argument("--fastq", "-fq", help='fastq file')
    parse.add_argument("--fasta", "-fa", help="output fasta, will erase if the file exists")
    args = parse.parse_args()
    assert args.fastq != args.fasta
    fastq_to_fasta(args.fastq, args.fasta)


if __name__ == "__main__":
    main()
