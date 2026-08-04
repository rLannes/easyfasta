"""This module manages all .fai index related functions: building, loading,
and querying by coordinate, including the Position helper type."""
from __future__ import annotations
from dataclasses import dataclass
from collections import deque

import typing
from collections.abc import Iterable
from typing_extensions import TypedDict, NotRequired
from typing import TextIO, Generator, Optional
from .common import fasta_iter, wrap_sequence, reverse_complement

import os
from pathlib import Path
import logging
log = logging.getLogger(__name__)

class line_length_checker():
    """
    Tracks the last few line lengths seen while scanning a fasta sequence
    to enforce the .fai constraint that all lines must be the same size
    except (optionally) the last one.
    """
    def __init__(self):
        """
        :return None: starts with an empty history
        """
        self.data = deque()

    def check(self, size):
        """
        record a new line length and check it is consistent with the
        previous ones (only the most recent line may differ in size).

        :param int size: the length of the line just read
        :return bool: False if a length mismatch is detected, True otherwise
        """
        self.data.appendleft(size)
        if len(self.data) > 3:
            self.data.pop()
        if len(self.data) < 3:
            return True
        if self.data[1] != self.data[2]:
            return False
        return True

    def reset(self):
        """
        :return None: clears the recorded history, for use at record boundaries
        """
        self.data.clear()

    def size(self):
        """
        :return int: the most recently recorded line length
        """
        return self.data.pop()



def make_fasta_multiline(fasta, out):
    """
    helper function take a fasta file and make it multiline.
    usefull because it garantee all lines are the same length (required by .fai format) except the last one

    :param str|Path fasta: the fasta file to reformat
    :param str|Path out: path to write the reformatted, multiline fasta to
    :return None: writes the reformatted fasta to out

    """

    with open(fasta) as fi, open(out, 'w') as fo:
        for p,s in fasta_iter(fi):
            fo.write(">{}\n{}\n".format(p, wrap_sequence(s, 100)))




def fasta_index_fai(fasta):
    """
    Build a .fai index file next to fasta by scanning it once.
    have undefined behaviour if empty line inside a sequence.

    :param str|Path fasta: the fasta file to index
    :raises AssertionError: if the fasta file is empty or its sequence lines
        are not all the same size (only the last line of a record may differ)
    :return None: writes a fasta.fai index file next to fasta
    """
    fasta = Path(fasta)
    index = fasta.with_name(fasta.name + ".fai")
    try:
        assert os.stat(fasta).st_size != 0
    except AssertionError:
        raise AssertionError("fasta file is empty ERROR")
    except:
        raise

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
@dataclass
class Position:
    """A genomic interval."""
    chr: str
    start: int
    end: int
    strand: str = "." 

    def __lt__(self, other):
        """
        :param Position other: another position on the same chr
        :return bool: True if this position's start is before other's
        """
        assert self["chr"] == other["chr"]
        return self["start"] < other["start"]

    def __gt__(self, other):
        """
        :param Position other: another position on the same chr
        :return bool: True if this position's start is after other's
        """
        assert self["chr"] == other["chr"]
        return self["start"] > other["start"]

    def __le__(self, other):
        """
        :param Position other: another position on the same chr
        :return bool: True if this position's start is before or equal to other's
        """
        assert self["chr"] == other["chr"]
        return self["start"] <= other["start"]

    def __ge__(self, other):
        """
        :param Position other: another position on the same chr
        :return bool: True if this position's start is after or equal to other's
        """
        assert self["chr"] == other["chr"]
        return self["start"] >= other["start"]

    @staticmethod
    def from_dict(dict):
        """
        :param dict dict: a dict-like object with keys "chr", "start", "end", and optionally "strand"
        :return Position: a new Position built from the given dict
        """
        return Position(dict["chr"], dict["start"], dict["end"], dict.get("strand", "."))


    def to_dict(self):
        return {"chr": self["chr"], "start": self["start"], "end": self["end"], "strand": self["strand"]}

    
    def __len__(self):
        """
        :return int: the number of fields in a Position (always 4)
        """
        return 4

    def size(self):
        """
        :return int: the length of the interval, end - start
        """
        return self["end"] - self["start"]

    def __eq__(self, other):
        """
        :param Position other: another position or dict-like object to compare against
        :return bool: True if chr, start, end, and strand all match. strand defaults
            to "." when absent from other, matching Position's own default
        """
        if self["chr"] != other["chr"]:
            return False
        if self["start"] != other["start"]:
            return False
        if self["end"] != other["end"]:
            return False
        if self["strand"] != other.get("strand", "."):
            return False
        return True


    def __getitem__(self, key):
        """
        :param str key: one of "chr", "start", "end", "strand"
        :return: the value stored under key
        :raises KeyError: if key is not a valid field name
        """
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(key)

    def get(self, key, default=None):
        """
        :param str key: one of "chr", "start", "end", "strand"
        :param default: value to return if key is not a valid field name
        :return: the value stored under key, or default
        """
        return getattr(self, key, default)

    def keys(self):
        """
        :return tuple[str]: the field names of a Position
        """
        return ("chr", "start", "end", "strand")

    def __contains__(self, key):
        """
        :param str key: a candidate field name
        :return bool: True if key is one of Position's field names
        """
        return key in self.keys()


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


    fasta = Path(fasta)
    index = fasta.with_name(fasta.name + ".fai")
    if not index.is_file():
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

    fasta = Path(fasta)
    index = fasta.with_name(fasta.name + ".fai")
    if not index.is_file():
        raise FileNotFoundError("cannot find the fai file, is your fasta indexed? ")

    if dico_index is None:
        dico_index = _get_dico_index(index)

    seqs = ""
    with open(fasta) as fi:
        fi.seek(0)
        for p in sorted(positions, key=lambda x: x["start"]):
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
    fasta = Path(fasta)
    index = fasta.with_name(fasta.name + ".fai")
    if not index.is_file():
        raise FileNotFoundError("cannot find the fai file, is your fasta indexed? ")

    if dico_index is None:
        dico_index = _get_dico_index(index)

    (start_bytes, end_bytes) = _get_bytes(name, dico_index, start, end)

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

