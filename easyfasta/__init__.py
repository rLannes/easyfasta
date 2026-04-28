from __future__ import annotations
"""
v1.1.0
Now support fai index due to popular demand.
Breaking changes from v1.0.14:
- build_index() now builds a .fai index file, use build_dico_index() for dictionary index
- get_sequence_index() has been renamed to get_sequence_dico_index()

v.1.2 
now come with fastq to fasta will soon implement fai for fastq
v.1.2.6
query function family extended with member accepting dict like object:
query_position
query_iter
query_splice

"""
from .common import(
fasta_iter,
complement,
reverse,
reverse_complement,
wrap_sequence,
)

from .easyfasta import(

    build_dico_index,
    get_sequence_dico_index,
    get_sequence_id,
    load_fasta,
    build_index,
    load_index,
    query,


)

from .fai_common import(
    query_position,
    query_iter,
    query_splice
)

__all__ = [
    "wrap_sequence",
    "build_index",
    "build_dico_index",
    "load_index",
    "fasta_iter",
    "query",
    "query_position",
    "query_iter",
    "query_splice",
    "get_sequence_dico_index",
    "get_sequence_id",
    "build_index",
    "fasta_iter",
    "load_fasta",
    "complement",
    "reverse",
    "reverse_complement"
]

__version__ = "1.0.5"

