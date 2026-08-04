"""
Low-level, format-agnostic building blocks: sequence complement/reverse
helpers, the FastaRecord container, and the fasta/fastq line iterators.
"""
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



class FastaRecord():
    """ 
    A FASTA record storing a description and sequence.
    Supports tuple-style unpacking (desc, seq = record) and
    dict-style access (record["seq"]) for backward compatibility.

    :param str description: the description line of the FASTA record
    :param str seq: the nucleotide or amino acid sequence

    :attribute str id: the first word of the description line
    :attribute str description: the full description line
    :attribute str seq: the full sequence
    :attribute int len: the length of the sequence

    :note: use str(record) or f"{record}" to get a valid FASTA string

    :example:

        >>> record = FastaRecord("gene1 some info", "ATCGATCG")
        >>> record.id
        'gene1'
        >>> record.description
        'gene1 some info'
        >>> record["seq"]
        'ATCGATCG'
        >>> desc, seq = record
    """
    def __init__(self, description, seq):
        """
        :param str description: the description line of the FASTA record
        :param str seq: the nucleotide or amino acid sequence
        """
        self._description = description.strip()
        self._seq = seq.strip()

    def __iter__(self):
        """
        :return Iterator[str, str]: (description, seq), for tuple-style unpacking
        """
        return iter((self._description, self._seq))

    def __getitem__(self, value):
        """
        dict-style or tuple-style access to the record's fields.

        :param int|str value: 0/1 for tuple-style access, or one of "seq", "id", "description"
        :return str: the requested field
        """
        if isinstance(value, int):
            assert 0 <= value <= 1
            return (self._description, self._seq)[value]

        assert value in ["seq", "id", "description"]
        return getattr(self, value)

    @property
    def len(self):
        """
        :return int: the length of the sequence
        """
        return len(self._seq)

    @property
    def id(self):
        """
        :return str: the first word of the description line
        """
        return self._description.split()[0]

    @property
    def description(self):
        """
        :return str: the full description line
        """
        return self._description

    @property
    def seq(self):
        """
        :return str: the full sequence
        """
        return self._seq

    def __eq__(self, other):
        """
        :param FastaRecord other: the record to compare against
        :return bool: True if description and seq are equal
        """
        return (self._description, self._seq) == (other._description, other._seq)

    def __str__(self):
        """
        :return str: a valid, wrapped FASTA-formatted string for this record
        """
        return f">{self._description}\n{wrap_sequence(self._seq)}"

    def __repr__(self):
        """
        :return str: a debug representation showing id and sequence length
        """
        return f"FastaRecord('{self.id}', len={self.len} bp)"


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



def fasta_iter(open_file: TextIO, position: bool=None) -> Generator[FastaRecord, None, None] |  Generator[tuple[FastaRecord, int], None, None]:
    """
    An Iterator over an opened fasta file.

    Note: I developed this while working on extremely large fasta files, 
    which make no sense to load into memory.

    Yields FastaRecord objects that support both tuple-style unpacking 
    and dict-style access.

    .. code-block:: python

        # tuple-style (backward compatible)
        with open(fasta_file) as fi:
            for description, sequence in fasta_iter(fi):
                sequence_id = description.split()[0]
                print(description, sequence_id, sequence)

        # FastaRecord-style
        with open(fasta_file) as fi:
            for record in fasta_iter(fi):
                print(record.id, record.seq, record.len)
                print(record["seq"])

        # with position tracking
        with open(fasta_file) as fi:
            for record, pos in fasta_iter(fi, position=True):
                print(record.id, pos)

    :note: you can use str(record) or f"{record}" to get a valid FASTA string

    :param TextIO open_file: an opened fasta file
    :param bool position: if True, also yield the byte offset of each record 
        as returned by tell()

    :yield FastaRecord: a fasta record (description, sequence)
    :yield tuple(FastaRecord, int): if position is True, a tuple of 
        (record, byte_offset)
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
                    yield FastaRecord(p, seq)
                else:
                    yield FastaRecord(p, seq), last_pos
                p, seq = "", ""
                last_pos = pos
                pos = open_file.tell()
                
            p = line[1:].strip()

        else:
            seq += line.strip()
            pos = open_file.tell()
        line = open_file.readline()

    if not position:
        yield FastaRecord(p, seq)
    else:
        yield FastaRecord(p, seq), last_pos


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


