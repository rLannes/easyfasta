from __future__ import annotations
"""This module manage all fai related functions"""
from collections import deque

import typing
from collections.abc import Iterable
from typing import TextIO, Generator
import common
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
        for p,s in common.fasta_iter(fi):
            fo.write(">{}\n{}\n".format(p, common.wrap_sequence(s, 100)))

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
            if new_line.startswith(">"):
                if name: 
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


def query(fasta, name, start, end, strand="+", dico_index=None):
    """if strand == "-" return the reverse complement"""
    index = fasta + ".fai"
    if not Path(index).is_file():
        raise AssertionError("cannot find the fai file, is your fasta indexed? ")
    
    if dico_index is None:
        dico_index = {}
        with open(index) as fi:
            for l in fi:
                spt = l.strip().split()
                dico_index[spt[0]] = spt[1:]
    if name not in dico_index:
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
    end_bytes =  offset + (line_number * line_bytes) + base_offset

    seq = ""
    with open(fasta) as fi:
        fi.seek(start_bytes)
        seq = fi.read(end_bytes - start_bytes)

    seq = "".join(seq.strip().split())
    if strand == "-":
        seq = common.reverse_complement(seq)
    return seq


# def query_open_fasta(fasta, name, start, end, strand="+", dico_index=None):
#     """if strand == "-" return the reverse complement"""
#     index = fasta + ".fai"
#     if not Path(index).is_file():
#         raise AssertionError("cannot find the fai file, is your fasta indexed? ")
    
#     if dico_index is None:
#         dico_index = {}
#         with open(index) as fi:
#             for l in fi:
#                 spt = l.strip().split()
#                 dico_index[spt[0]] = spt[1:]
#     if name not in dico_index:
#         log.error("{} not in index".format(name))
#         raise AssertionError
#     length,offset,linebases,line_bytes = [int(x) for x in dico_index[name]]

#     if start > length or end > length or start < 0 or end < start:
#         log.error("coordinate error")
#         raise AssertionError


#     line_number = start // linebases
#     base_offset = start % linebases
#     start_bytes =  offset + (line_number * line_bytes) + base_offset
#     line_number = end // linebases
#     base_offset = end % linebases
#     end_bytes =  offset + (line_number * line_bytes) + base_offset

#     seq = ""
#     fasta.seek(0)
#     fi.seek(start_bytes)
#     seq = fi.read(end_bytes - start_bytes)

#     seq = "".join(seq.strip().split())
#     if strand == "-":
#         seq = easyfata.reverse_complement(seq)
#     return seq

