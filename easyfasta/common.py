from __future__ import annotations
import typing
from collections.abc import Iterable
from typing import TextIO, Generator

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


def fastq_iter(open_file: TextIO, position: bool=None)-> Generator[tuple[str, str, str], None, None] |  Generator[tuple[str, str, str, int], None, None]:
    """
    An Iterator over an opened FASTQ file.

    Note: Developed for large FASTQ files that should not be loaded into memory.
    Validates record structure by checking the '+' separator and that len(sequence) == len(quality).

    .. code-block:: python

        with open(fastq_file) as fi:
            for identifier, sequence, quality in fastq_iter(fi):
                print(identifier, sequence, quality)

    :param TextIO open_file: an opened FASTQ file
    :param bool position: if true return the byte offset of the record start as reported by tell. The signature becomes Generator((str, str, str, int))
    :return Generator((str, str, str)): Iterable(identifier, sequence, quality)
    """
    
    pos = 0
    
    id_, seq, qual = "", "", ""
    open_file.seek(0)

    line = open_file.readline()
    while line:

        if line.startswith('@'):
            id_ = line.strip()[1:]
            seq = open_file.readline().strip()
            if "+" != open_file.readline().strip():
                raise AssertionError("record id_: {}, seq:{}, at pos {} is broken".format(id_, seq, pos ))
            qual =  open_file.readline().strip()
            if len(qual) != len(seq):
                raise AssertionError("record id_: {}, seq:{}, qual: {}, at pos {}. is broken len(qual) != len(seq)".format(id_, seq,qual, pos,  ))
                                 
            if not position:
                 yield id_, seq, qual
            else:
                yield id_, seq, qual, pos
        pos = open_file.tell()
        line = open_file.readline()



def fasta_iter(open_file: TextIO, position: bool=None) -> Generator[tuple[str, str], None, None] |  Generator[tuple[str, str, int], None, None]:
    """
    An Iterator over an opened fasta file.

    Note: I developed this while working on extremely large fasta file, which make no sense to load into memory.

    .. code-block:: python

        with open(fasta_file) as fi:
            for identifier_line, sequence in fasta_iter(fi):
                sequence_id = identifier_line.split()[0]
                print(identifier_line, sequence_id, sequence)


   
    :param TextIO  open_file: an opened fasta file
    :param bool position: if true return the start of the sequence (including the identifier line) returned by tell. and the signature become Generator((str, str, int)) 
    :return Generator((str, str)): Iterable(prompt, sequence)
    """

    pos = 0
    last_pos = 0
    
    p, seq = "", ""
    open_file.seek(0)

    line = open_file.readline()
    while line:

        if line.startswith('>'):

            if seq:

                if not position:
                    yield p, seq
                else:
                    yield p, seq, last_pos
                p, seq = "", ""
                last_pos = pos
                pos = open_file.tell()
                
            p = line[1:].strip()

        else:
            seq += line.strip()
            pos = open_file.tell()
        line = open_file.readline()

    if not position:
        yield p, seq
    else:
        yield p, seq, last_pos

def wrap_sequence(sequence: str, chunk_size: int=80) -> str:
    """
    chunk a string in multiple lines by adding a new line every chunk size
    useful to write multiline fasta.

    :param str sequence: the string to make multiline
    :param int chunk_size: the size of the line.

    :return str: 
    
    """
    new_seq = ""
    cpt = 0
    while cpt <= len(sequence):
        new_seq += sequence[cpt: cpt + chunk_size] + "\n"
        cpt += chunk_size
    return new_seq.strip()