import typing
from pathlib import Path
from collections.abc import Iterable
from typing import TextIO, Generator


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


def build_index(fasta_file: str|Path) -> dict[str, int]:
    """
    build an index from a fasta file, dict sequence identifier -> position

    :param str|Path fasta_file: the fasta file to build index from
    :return  dict[str, int]: index dictionary  identifier -> position
    """
    index = {}
    with open(fasta_file) as fi:
        for p, s, i in fasta_iter(fi, position=True):
            index[p.split()[0]] = i
    return index

def get_sequence_index(fasta_file: str|Path, identifiers:Iterable[str], index_dict:dict[str, int], ignore_unfound: bool = True) -> list[tuple[str, str]]:
    """
    uses index to get sequence from a file faster than just parsing through the file. you need to generate an index first (you can use build_index)
    will not raise an error if any identifier in identifiers are not in the dict. you can turn off this by setting ignore_unfound to True


    .. code-block:: python
        index = build_index(fasta_file)
        # you can save/load the index using pickle 
        pickle.dump(index, filename)
        index = pickle.load(filename)
        # this can save large amount of time on large file
        sequences = get_sequence_index(fasta_file, identifiers, index)

    :param str|Path  fasta_file: a fasta file
    :param Iterable identifier: an iterable with id to recover sequence from
    :param dict[str, int] index_dict: a dictionary associating identifier to a position in file, you can make one from build_index
    :param bool ignore_unfound: defualt False.
    :return [(str, str)]: [(identifier, sequence)] for each sequences with identifier present in identifier

    """
    res = []
    with open(fasta_file) as open_file:

        for id_ in identifiers:

            offset = index_dict.get(id_)
            if offset is None and ignore_unfound:
                continue
            elif offset is None:
                print("id: {} is not in index".format(id_))
                index_dict[id_]  #raise eror
            open_file.seek(offset, 0)

            open_file.readline()
            
            sequence = ""
            line = open_file.readline().strip()

            while not line.startswith(">") and line:
                sequence += line.strip()
                line = open_file.readline()
            res.append((id_, sequence))
    return res


def get_sequence_id(fasta_file: str|Path, identifiers: Iterable[str], identifier_only: bool=True) -> list[tuple[str, str]]:
    """
    return sequence in identifiers from the fasta_file. !! will NOT throw a warning/error if a sequence is not found in the fasta!!

    :param str|Path  fasta_file: a fasta file
    :param Iterable identifier: an iterable with id to recover sequence from
    :param bool identifier_only: fasta are composed of identifier and metadata, by default only use the identifier part of the fasta line set to false to use the full line.
    :return [(str, str)]: [(identifier, sequence)] for each sequences with identifier present in identifier

    """
    res = []

    with open(fasta_file) as open_file:
        for p, s in fasta_iter(open_file=open_file):
            if identifier_only:
                p = p.split()[0]
            if p in identifiers:
                res.append((p, s))
    
    return res

# TODO change to fasta iterator
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


def load_fasta(fasta) -> dict[str, str]:
    """
    return dictionary association sequence identifier to its sequence from a fasta file 
    
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
