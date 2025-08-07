# FASTA Utils

A lightweight functional Python library for efficient FASTA file parsing and DNA sequence manipulation. No OOP bloat, data only data.

## Features

- **Memory-efficient parsing**: Stream through large FASTA files without loading everything into memory
- **Random access**: Jump directly to specific sequences with position tracking
- **Sequence extraction**: Filter sequences by identifiers
- **DNA manipulation**: Complete IUPAC-compliant complement and reverse complement operations
- **Formatting**: Convert sequences to multi-line FASTA format
- **Does not validate input**: user are responsible to provide correctly formatted file.
## Installation

Simply copy the module to your project or install via pip (if published).

## Quick Start

```python
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

# Extract specific sequences using indexes
index = build_index('sequences.fasta')
# using pickle you can save and load the index
#pickle.dump(index, "save_index_file.pkl")
#index = pickle.load("save_index_file.pkl")
target_ids = ['seq1', 'seq2', 'seq3']
found = get_sequence_index(f, target_ids, index, ignore_unfound=True)
for header, seq in found:
    print(f"Found: {header}")

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

#### `build_index(fasta_file: str|Path) -> dict[str, int]`

Build a fasta index as a dictionary


```python
index = build_index(fasta_file)
```

#### `get_sequence_index(fasta_file: str|Path, identifiers:Iterable[str], index_dict:dict[str, int], ignore_unfound: bool = False) -> list[tuple[str, str]]`

use index to retrieve sequence (faster)


```python
index = build_index(fasta_file)
wanted = ['seq1', 'seq2']
get_sequence_index(fasta_file, wanted, index)
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
   fo.write(">{}\n{}\n".format('seq_id',  wrap_sequence("ATCGATCGATCG" * 10, 80)))
```

## Design Philosophy

This library prioritizes:

- **Memory efficiency**: Built for large genomic files that don't fit in RAM
- **Simplicity**: Clean, predictable API with minimal dependencies. Not OOP bloat, only data.
- **Performance**: Stream-based processing with O(1) memory usage for parsing
- **Standards compliance**: Full IUPAC nucleotide code support

## Use Cases

- Processing large genome assemblies
- Extracting specific genes or regions
- Converting between sequence formats
- Quality control pipelines
- Bioinformatics workflows requiring memory efficiency

## Requirements

- Python 3.8+
- No external dependencies

## License
MIT

## Contributing
Feel free to ask for new features. I published it as lightweight because those are the feature I use the most and wanted to start with a solid fondation.

I used this library for years, and it has been extensively tested. As such I will only adress issue that come with a minimal reproducible problem.
