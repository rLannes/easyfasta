import typing
from pathlib import Path
from collections.abc import Iterable
from typing import TextIO, Generator


def chunk_seq(sequence: str, chunk_size: int=80) -> str:
    """
    chunk a string in mulitple line by adding a new line every chunck size
    usefull to write mulit line fasta.

    :param str sequence: the string to make multiline
    :param int chunksize: the size of the line.
    
    """
    new_seq = ""
    cpt = 0
    while cpt <= len(sequence):
        new_seq += sequence[cpt: cpt + chunk_size] + "\n"
        cpt += chunk_size
    return new_seq.strip()


def extract_sequences(open_file: TextIO, identifiers: Iterable[str], identifier_only: bool=True) -> list[tuple[str, str]]:
    """
    return sequence in identifiers from the opened fasta_file open_file. !! will NOT throw a warning/error if a sequence is not found in the fasta!!

    :param TextIO  open_file: an opened fasta file
    :param Iterable identifier: an iterable with id to recover sequence from
    :param bool identifier_only: fasta are composed of identifier and metadata, by default only use the identifier part of the fasta line set to false to use the full line.
    :return [(str, str)]: [(identifier, sequence)] for each sequences with identifier present in identifier

    """
    res = []
    open_file.seek(0)
    for p, s in fasta_iter(open_file=open_file):
        if identifier_only:
            p = p.split()[0]
        if p in identifiers:
            res.append((p, s))
    
    return res

# TODO change to fasta iterator
def fasta_iter(open_file: TextIO, position: bool=None) -> Generator[tuple[str, str], None, None]:
    """
    An Iterator over an opened fasta file.

    Note: I developed this while working on extremly large fasta file, which make no sense to load into memory.

    .. code-block:: pyhton

        with open(fasta_file) as fi:
            for identifier_line, sequence in fasta_iter(fi):
                sequence_id = identifier_line.split()[0]
                print(identifier_line, sequence_id, sequence)


   
    :param TextIO  open_file: an opened fasta file
    :param bool position: if true return the start of the sequence (including the identifier line) returned by tell. and the signature become Generator((str, str, int)) 
    :return Generator((str, str)): Iterable(prompt, sequence)
    """

    pos = 0
    p, seq = "", ""
    open_file.seek(0)

    line = open_file.readline()
    while line:

        if line.startswith('>'):

            if seq:

                if not position:
                    yield p, seq
                else:
                    yield p, seq, pos
                p, seq = "", ""

            p = line[1:].strip()

        else:
            seq += line.strip()
            pos = open_file.tell()
        line = open_file.readline()

    if not position:
        yield p, seq
    else:
        yield p, seq, pos


def load_fasta(fasta) -> dict[str: str]:
    """
    return dictionnary association sequence identifier to its sequence from a fasta file 
    
    :param str|Path fasta: a fasta file
    :return dict[str: str]: identifier => sequence
    
    """
    result = {}
    with open(fasta) as fi:
        for p, s in fasta_iter(fi):
            result[p.split()[0]] = s
    return result

""" constant IUPAC complete DNA complement, case insensitive"""
DNA_COMPLEMENT={"A":"T", "T":"A", "C":"G", "G":"C", "N":"N", "S":"S",\
                "W":"W", "R": "Y", "Y": "R", "K": "M", "M": "K", "B":"V", "V":"B", "D":"H", "H":"D",\
                "s":"s", "w":"w", "r": "y", "y": "r", "k": "m", "m": "k", "b":"v", "v":"b", "d":"h", "h":"d",\
                "a":"t", "t":"a", "c":"g", "g":"c", "n":"n"}

def complement(seq: str) -> str:
    """
    case insensitive IUPAC complete complement of a DNA sequence

    :param str seq: a DNA sequence
    :return str: Complement sequence

    """
    return "".join([DNA_COMPLEMENT[x] for x in seq ])

def reverse(seq: str) -> str:
    """
    reverse the given sequence

    :param str seq: a DNA sequence
    :return str: reverse sequence

    """
    return "".join([x for x in seq[::-1]])

def reverse_complement(seq: str) -> str:
    """
    case insensitive IUPAC complete reverse complement of a DNA sequence

    :param str seq: a DNA sequence
    :return str: reverse complemented sequence

    """
    return "".join([DNA_COMPLEMENT[x] for x in seq[::-1]])



