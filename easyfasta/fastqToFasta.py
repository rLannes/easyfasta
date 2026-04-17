import argparse
import typing
from collections.abc import Iterable
from typing import TextIO, Generator
from .common import fastq_iter


def fastq_to_fasta(fastq: str, fasta: str):
    with open(fastq) as fi, open(fasta, "w") as fo:
        for record in fastq_iter(fi):
            fo.write(">{}\n{}\n".format(record[0], record[1]))

def main():
    parse = argparse.ArgumentParser(description="make a fasta from a fastq")
    parse.add_argument("--fastq", "-fq", help='fastq file')
    parse.add_argument("--fasta", "-fa", help="output fasta, will erase if the file exists")
    args = parse.parse_args()
    assert args.fastq != args.fasta
    fastq_to_fasta(args.fastq, args.fasta)
    

if __name__ == "__main__":
    main()
