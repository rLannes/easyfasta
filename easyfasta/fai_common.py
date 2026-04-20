from __future__ import annotations
"""This module manage all fai related functions"""
from collections import deque

import typing
from collections.abc import Iterable
from typing_extensions import TypedDict, NotRequired
from typing import TextIO, Generator, Optional
from .common import fasta_iter, wrap_sequence, reverse_complement


from pathlib import Path
import logging
log = logging.getLogger(__name__)

class line_length_checker():
    def __init__(self):
        self.data = deque()

    def check(self, size):
        self.data.appendleft(size)
        if len(self.data) > 3:
            self.data.pop()
        if len(self.data) < 3:
            return True
        if self.data[1] != self.data[2]:
            return False
        return True

    def reset(self):
        self.data.clear()

    def size(self):
        return self.data.pop()



def make_fasta_multiline(fasta, out):
    """
    helper function take a fasta file and make it multiline.
    usefull because it garantee all lines are the same length (required by .fai format) except the last one

    """

    with open(fasta) as fi, open(out, 'w') as fo:
        for p,s in fasta_iter(fi):
            fo.write(">{}\n{}\n".format(p, wrap_sequence(s, 100)))

def fasta_index_fai(fasta):  

    index = fasta + ".fai" 
    name = None
    offset = None
    seq_size = None
    line_base_checker = line_length_checker()
    line_off_checker = line_length_checker()
    with open(fasta) as fi, open(index, "w") as fo:
        prev_line_cursor = 0
        new_line = fi.readline()

        while new_line:
            if not new_line.strip(): # empty line
                new_line = fi.readline()
                continue
                
            if new_line.startswith(">"):
                if name is not None: 
                    fo.write("{}\t{}\t{}\t{}\t{}\n".format(name, str(seq_size), str(offset), str(line_base_checker.size()), str(line_off_checker.size())))
                line_base_checker.reset()
                line_off_checker.reset()
                name = new_line.split()[0][1:]
                offset = fi.tell()
                seq_size = 0

            else:
                base_l = len(new_line.strip())
                seq_size += base_l
                if not line_base_checker.check(base_l):
                    raise AssertionError("Sequence line must be the same size, only the last line can differ in size. You can use make_fasta_multiline to proof yout fasta file")
                offset_line = fi.tell() - prev_line_cursor 
                if not line_off_checker.check(offset_line):
                    raise AssertionError("Sequence line must be the same size, this error was likely trigger by a whitespace of invisible characters" \
                                        " only the last line can differ in size")

        
            prev_line_cursor = fi.tell()
            new_line = fi.readline()

        fo.write("{}\t{}\t{}\t{}\t{}\n".format(name, str(seq_size), str(offset), str(line_base_checker.size()), str(line_off_checker.size())))



# class for typeHint
class Position(TypedDict):
    """A genomic interval."""
    chr: str
    start: int
    end: int
    strand: NotRequired[str]

def query_position(fasta: str, position: Position,
                   dico_index: Optional[dict] = None) -> str:
    """
    Convenience wrapper around query() that accepts a dict-like genomic position.
    position = TypedDict("Position", {"chr": str, "start": int, "end": int, "strand": NotRequired[str]})
    Coordinates are zero-based (start inclusive, end exclusive).


    .. code-block:: python

        pos = {"chr": "chr1", "start": 100, "end": 200, "strand": "-"}
        seq = query_position("genome.fa", pos)

    :param str fasta: path to an indexed fasta file (.fai must exist)
    :param position: any object supporting dict-style access with keys
        "chr" (str), "start" (int), "end" (int), and optionally
        "strand" (str, defaults to "+")
    :param dict dico_index: pre-loaded index dict, passed through to query()
    :return str: the extracted sequence (reverse-complemented if strand is "-")
    :seealso: query()
    """
    return query(fasta, position["chr"], position["start"], 
                 position["end"], position.get("strand", "+"), dico_index) 


def query_iter(fasta: str, positions: list[Position],
               dico_index: Optional[dict] = None) -> list[str]:
    """
    Query multiple independent regions from an indexed fasta file.
    position = TypedDict("Position", {"chr": str, "start": int, "end": int, "strand": NotRequired[str]})

    Returns one sequence per position. The file is opened once and seeked
    for each region. Reverse complement is applied per-position when
    ``strand`` is ``"-"``.

    Coordinates are zero-based (start inclusive, end exclusive).

    .. code-block:: python

        positions = [
            {"chr": "chr1", "start": 100, "end": 200},
            {"chr": "chr2", "start": 500, "end": 600, "strand": "-"},
        ]
        seqs = query_iter("genome.fa", positions)

    :param str fasta: path to an indexed fasta file (.fai must exist)
    :param list positions: list of dict-like objects with keys
        "chr" (str), "start" (int), "end" (int), and optionally
        "strand" (str, defaults to "+")
    :param dict dico_index: pre-loaded index dict. If None, read from disk.
    :raises AssertionError: if the .fai file is missing or coordinates
        are invalid
    :return list[str]: one sequence per position
    :see also: query_splice() for concatenated exon queries
    """


    index = fasta + ".fai"
    if not Path(index).is_file():
        raise FileNotFoundError("cannot find the fai file, is your fasta indexed? ")
    
    if dico_index is None:
        dico_index = _get_dico_index(index)

    seqs = []
    with open(fasta) as fi:
        fi.seek(0)
        for p in positions:
            (start_bytes, end_bytes) = _get_bytes( p["chr"], dico_index,
                                                   p["start"], p["end"])

            fi.seek(start_bytes)
            this_seq = fi.read(end_bytes - start_bytes)
            this_seq = "".join(this_seq.strip().split())
            if p.get("strand", "+") == "-":
                this_seq = reverse_complement(this_seq)
            seqs.append(this_seq)
    return seqs


def query_splice(fasta: str, positions: list[Position],
                 dico_index: Optional[dict] = None) -> str:
    """
    Query and concatenate multiple regions (e.g. exons) into a single sequence.
    position = TypedDict("Position", {"chr": str, "start": int, "end": int, "strand": NotRequired[str]})

    All positions are read in order and joined. If the strand of the
    **first** position is ``"-"``, the concatenated result is
    reverse-complemented. Mixed strands (trans-splicing) are not supported
    and will produce undefined behaviour.

    Coordinates are zero-based (start inclusive, end exclusive).

    .. code-block:: python

        exons = [
            {"chr": "chr1", "start": 1000, "end": 1200, "strand": "+"},
            {"chr": "chr1", "start": 1500, "end": 1700, "strand": "+"},
        ]
        transcript_seq = query_splice("genome.fa", exons)

    :param str fasta: path to an indexed fasta file (.fai must exist)
    :param list positions: list of dict-like objects with keys
        "chr" (str), "start" (int), "end" (int), and optionally
        "strand" (str, defaults to "+"). All entries should share the
        same strand.
    :param dict dico_index: pre-loaded index dict. If None, read from disk.
    :raises AssertionError: if the .fai file is missing or coordinates
        are invalid
    :return str: the concatenated sequence, reverse-complemented if
        the first position's strand is "-"
    :see also: query_iter() for independent per-region queries
    """

    index = fasta + ".fai"
    if not Path(index).is_file():
        raise FileNotFoundError("cannot find the fai file, is your fasta indexed? ")
    
    if dico_index is None:
        dico_index = _get_dico_index(index)

    seqs = ""
    with open(fasta) as fi:
        fi.seek(0)
        for p in positions:
            (start_bytes, end_bytes) = _get_bytes( p["chr"], dico_index,
                                                   p["start"], p["end"])

            fi.seek(start_bytes)
            this_seq = fi.read(end_bytes - start_bytes)
            seqs +=  "".join(this_seq.strip().split())
        
    if positions[0].get("strand", "+") == "-":
        return reverse_complement(seqs)
    return seqs


def query(fasta: str, name: str, start: int, end: int, strand="+", dico_index=None) -> str:
    """
    Query a sequence from an indexed fasta file using the .fai index.
    Uses byte offset seeking, so memory usage is independent of file size.
    note: position are zero based

    If strand is "-", returns the reverse complement.

    .. code-block:: python

        seq = query("genome.fa", "chr1", 100, 200)
        seq = query("genome.fa", "chr1", 100, 200, strand="-")

        # reuse the index for multiple queries
        idx = None
        for name, start, end in positions:
            seq = query("genome.fa", name, start, end, dico_index=idx)

    :param str fasta: path to an indexed fasta file (.fai must exist)
    :param str name: sequence/chromosome name as found in the .fai index
    :param int start: start position (0-based, inclusive)
    :param int end: end position (0-based, exclusive)
    :param str strand: "+" (default) or "-" for reverse complement
    :param dict dico_index: pre-loaded index dict to avoid re-reading the 
        .fai file on repeated queries. If None, the index is read from disk.

    :raises AssertionError: if the .fai file is missing, the name is not 
        in the index, or coordinates are out of bounds

    :return str: the extracted sequence
    """
    index = fasta + ".fai"
    if not Path(index).is_file():
        raise FileNotFoundError("cannot find the fai file, is your fasta indexed? ")
    
    if dico_index is None:
        dico_index = _get_dico_index(index)

    (start_bytes, end_bytes) = _get_bytes(name, dico_index, start, end)

    """if name not in dico_index:
        log.error("{} not in index".format(name))
        raise AssertionError
    length,offset,linebases,line_bytes = [int(x) for x in dico_index[name]]

    if start > length or end > length or start < 0 or end < start:
        log.error("coordinate error")
        raise AssertionError


    line_number = start // linebases
    base_offset = start % linebases
    start_bytes =  offset + (line_number * line_bytes) + base_offset
    line_number = end // linebases
    base_offset = end % linebases
    end_bytes =  offset + (line_number * line_bytes) + base_offset"""

    seq = ""
    with open(fasta) as fi:
        fi.seek(start_bytes)
        seq = fi.read(end_bytes - start_bytes)

    seq = "".join(seq.strip().split())
    if strand == "-":
        seq = reverse_complement(seq)
    return seq


def _get_bytes(name, dico_index, start, end):
    """
    Convert genomic coordinates to byte offsets in the fasta file.

    :param str name: sequence/chromosome name
    :param dict dico_index: index dict as returned by _get_dico_index()
    :param int start: start position (0-based, inclusive)
    :param int end: end position (0-based, exclusive)
    :raises AssertionError: if *name* is not in the index or coordinates
        are out of bounds
    :return tuple[int, int]: (start_bytes, end_bytes) file offsets
    """
    if name not in dico_index:
        log.error("{} not in index".format(name))
        raise ValueError("{} not in index".format(name))
    length,offset,linebases,line_bytes = [int(x) for x in dico_index[name]]

    if start > length or end > length or start < 0 or end < start:
        log.error("coordinate error")
        raise ValueError("coordinate error")


    line_number = start // linebases
    base_offset = start % linebases
    start_bytes =  offset + (line_number * line_bytes) + base_offset
    line_number = end // linebases
    base_offset = end % linebases
    end_bytes =  offset + (line_number * line_bytes) + base_offset

    return (start_bytes, end_bytes)

def _get_dico_index(index):
    """
    Parse a ``.fai`` index file into a dict.

    :param str index: path to the .fai file
    :raises AssertionError: if the file does not exist
    :return dict: mapping of sequence name → list of index fields
        [length, offset, linebases, line_bytes]
    """
    if not Path(index).is_file():
        raise FileNotFoundError("cannot find the fai file, is your fasta indexed? ")

    dico_index = {}
    with open(index) as fi:
        for l in fi:
            spt = l.strip().split()
            dico_index[spt[0]] = spt[1:]

    return dico_index

