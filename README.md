# Easy Fasta

A lightweight functional Python library for efficient FASTA file parsing and DNA sequence manipulation. No OOP bloat, only data.

## Features

- **Memory-efficient parsing**: Stream through large FASTA files without loading everything into memory
- **Random access**: Jump directly to specific sequences with position tracking
- **FAI indexing**: Build and query standard `.fai` index files for fast random access
- **Sequence extraction**: Filter sequences by identifiers
- **DNA manipulation**: Complete IUPAC-compliant complement and reverse complement operations
- **Formatting**: Convert sequences to multi-line FASTA format
- **Does not validate input**: users are responsible to provide correctly formatted files.

## new in 1.2.6

query function accepting dict like object:
query_position
query_iter
query_splice


## Installation
python 3.8+
```
> pip install easyfasta
```
or simply copy the module to your project

## Quick Start

```python
from easyfasta import *

# Parse FASTA file sequence by sequence (memory efficient)
with open('sequences.fasta') as f:
    for header, sequence in fasta_iter(f):
        print(f">{header}")
        print(sequence[:50])  # First 50 bases

# Load entire FASTA into dictionary
sequences = load_fasta('sequences.fasta')
print(sequences['sequence_id'])

# Extract specific sequences
target_ids = ['seq1', 'seq2', 'seq3']
found = get_sequence_id('sequences.fasta', target_ids)
for header, seq in found:
    print(f"Found: {header}")

# Extract specific sequences using a dictionary index
index = build_dico_index('sequences.fasta')

target_ids = ['seq1', 'seq2', 'seq3']
found = get_sequence_dico_index('sequences.fasta', target_ids, index, ignore_unfound=True)
for header, seq in found:
    print(f"Found: {header}")

# FAI index for fast random access
build_index('sequences.fasta')  # creates sequences.fasta.fai
index = load_index('sequences.fasta')  # load into memory for repeated queries
seq = query('sequences.fasta', 'seq1', 0, 100, strand='+', dico_index=index)
# or 
seq = query_position('sequences.fasta', {"chr":'seq1', "start": 0, "end": 100, "strand"='+'}, dico_index=index)

# DNA manipulation
dna = "ATCGGTAA"
print(complement(dna))           # TAGCCATT
print(reverse_complement(dna))   # TTACCGAT
```

## API Reference

### Parsing Functions

#### `fasta_iter(open_file: TextIO) -> Generator[tuple[str, str], None, None]`

Memory-efficient iterator over FASTA sequences.

```python
with open('large_file.fasta') as f:
    for header, sequence in fasta_iter(f):
        # Process one sequence at a time
        process_sequence(header, sequence)
```

#### `load_fasta(fasta_path: str|Path) -> dict[str, str]`

Load entire FASTA file into a dictionary mapping sequence IDs to sequences.

```python
sequences = load_fasta('sequences.fasta')
my_sequence = sequences['sequence_id']
```

#### `get_sequence_id(fasta_file: str|Path, identifiers: Iterable[str], identifier_only: bool = True) -> list[tuple[str, str]]`

Extract sequences matching specific identifiers.

- `identifier_only`: If True, match only the first part of headers (before whitespace)

```python
wanted = ['seq1', 'seq2']
results = get_sequence_id('sequences.fasta', wanted)
```

### Dictionary Index Functions

#### `build_dico_index(fasta_file: str|Path) -> dict[str, int]`

Build an in-memory index as a dictionary mapping sequence identifiers to their byte position in the file.

```python
index = build_dico_index('sequences.fasta')
```

#### `get_sequence_dico_index(fasta_file: str|Path, identifiers: Iterable[str], index_dict: dict[str, int], ignore_unfound: bool = True) -> list[tuple[str, str]]`

Use a dictionary index to retrieve sequences faster than parsing through the file.

```python
index = build_dico_index('sequences.fasta')
wanted = ['seq1', 'seq2']
results = get_sequence_dico_index('sequences.fasta', wanted, index)
```

### FAI Index Functions

#### `build_index(fasta: str|Path) -> None`

Build a standard `.fai` index file next to the fasta file. Required before using `load_index` or `query`.

```python
build_index('sequences.fasta')  # creates sequences.fasta.fai
```

#### `load_index(fasta: str|Path) -> dict[str, list]`

Load a `.fai` index file into memory for repeated queries.

```python
index = load_index('sequences.fasta')
```

#### `query(fasta: str|Path, name: str, start: int, end: int, strand: str = "+", dico_index: dict = None) -> str`

Query a fasta file for a sequence by name and coordinates using the FAI index. Returns the reverse complement if strand is `"-"`.

```python
build_index('sequences.fasta')
index = load_index('sequences.fasta')
seq = query('sequences.fasta', 'chr1', 1000, 2000, strand='+', dico_index=index)
```

### Sequence Manipulation

#### `complement(seq: str) -> str`
Return the complement of a DNA sequence (A↔T, C↔G, supports all IUPAC codes).

#### `reverse(seq: str) -> str`
Return the reverse of a sequence.

#### `reverse_complement(seq: str) -> str`
Return the reverse complement of a DNA sequence.

#### `wrap_sequence(sequence: str, chunk_size: int = 80) -> str`
Format sequence with line breaks every `chunk_size` characters (standard multiline FASTA format).

```python
formatted = wrap_sequence("ATCGATCGATCG" * 10, 60)
print(formatted)  # 60 characters per line
# write to a file
with open(out_file, 'w') as fo:
   fo.write(">{}\n{}\n".format('seq_id', wrap_sequence("ATCGATCGATCG" * 10, 80)))
```

## Migration Guide: 1.0.14 → 1.1.0

Version 1.1.0 introduces FAI index support and contains **breaking changes**.

### Breaking Changes

| 1.0.14 | 1.1.0 | Notes |
|--------|-------|-------|
| `build_index()` | `build_dico_index()` | `build_index()` now builds a `.fai` file, not a dictionary |
| `get_sequence_index()` | `get_sequence_dico_index()` | straight rename |

### New in 1.1.0

- `build_index()` — builds a standard `.fai` index file
- `load_index()` — loads a `.fai` index into memory
- `query()` — fast random access to any sequence region by coordinates

### What you need to change

```python
# 1.0.14
index = build_index('sequences.fasta')
results = get_sequence_index('sequences.fasta', ids, index)

# 1.1.0
index = build_dico_index('sequences.fasta')
results = get_sequence_dico_index('sequences.fasta', ids, index)
```

> ⚠️ **Important**: `build_index()` no longer returns a dictionary. Calling it expecting a dictionary index will silently produce wrong results. Use `build_dico_index()` instead.

## Design Philosophy

This library prioritizes:

- **Memory efficiency**: Built for large genomic files that don't fit in RAM
- **Simplicity**: Clean, predictable API with minimal dependencies. Not OOP bloat, only data.
- **Performance**: Stream-based processing with O(1) memory usage for parsing
- **Standards compliance**: Full IUPAC nucleotide code support

## Use Cases

- Processing large fasta files (metagenome)
- Common DNA sequence manipulation
- Common operations on fasta including parsing, indexing and sequence retrieval
- Bioinformatics workflows requiring memory efficiency

## Requirements

- Python 3.8+
- No external dependencies

## License
MIT

## Contributing
Feel free to ask for new features. I published it as lightweight because those are the features I use the most and wanted to start with a solid foundation.

I used this library for years, and it has been extensively tested. As such I will only address issues that come with a minimal reproducible problem.