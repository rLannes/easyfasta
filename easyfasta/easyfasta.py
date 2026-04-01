from __future__ import annotations
import typing
from  easyfasta import fai_common, common
from pathlib import Path
from collections.abc import Iterable
from typing import TextIO, Generator
import logging
log = logging.getLogger(__name__)




def build_index(fasta):
    """
    Build a .fai index file for a fasta file.
    :param str|Path fasta: the fasta file to index
    :return None: creates a .fai file next to the fasta file
    """
    fai_common.fasta_index_fai(fasta)
    

def query(fasta, name, start, end, strand="+", dico_index=None):
    """
    Query a fasta file for a sequence by name and coordinates.
    :param str|Path fasta: the fasta file to query
    :param str name: the sequence identifier to query
    :param int start: the start position (0-based)
    :param int end: the end position (0-based, exclusive)
    :param str strand: the strand to query, "+" or "-", default "+"
    :param dict|None dico_index: a preloaded index dictionary, if None the index is loaded from disk
    :return str: the queried sequence, reverse complemented if strand is "-"
    """

    return fai_common.query(fasta, name, start, end, strand=strand, dico_index=dico_index)

def load_index(fasta):
    """
    Load a .fai index file into a dictionary for repeated queries.
    :param str|Path fasta: the fasta file whose index to load
    :return dict[str, list]: index dictionary identifier -> [length, offset, linebases, line_bytes]
    """
    index = fasta + ".fai"
    if not Path(index).is_file():
        raise AssertionError("cannot find the fai file, is your fasta indexed? ")
    
    if dico_index is None:
        dico_index = {}
        with open(index) as fi:
            for l in fi:
                spt = l.strip().split()
                dico_index[spt[0]] = spt[1:]
    return dico_index


def build_dico_index(fasta_file: str|Path) -> dict[str, int]:
    """
    build an index from a fasta file, dict sequence identifier -> position

    :param str|Path fasta_file: the fasta file to build index from
    :return  dict[str, int]: index dictionary  identifier -> position
    """
    index = {}
    with open(fasta_file) as fi:
        for p, s, i in common.fasta_iter(fi, position=True):
            index[p.split()[0]] = i
    return index

def get_sequence_dico_index(fasta_file: str|Path, identifiers:Iterable[str], index_dict:dict[str, int], ignore_unfound: bool = True) -> list[tuple[str, str]]:
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
        for p, s in common.fasta_iter(open_file=open_file):
            if identifier_only:
                p = p.split()[0]
            if p in identifiers:
                res.append((p, s))
    
    return res


def load_fasta(fasta) -> dict[str, str]:
    """
    return dictionary association sequence identifier to its sequence from a fasta file 
    
    :param str|Path fasta: a fasta file
    :return dict[str: str]: identifier => sequence
    
    """
    result = {}
    with open(fasta) as fi:
        for p, s in common.fasta_iter(fi):
            result[p.split()[0]] = s
    return result



