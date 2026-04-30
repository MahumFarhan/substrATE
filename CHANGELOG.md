# Changelog

All notable changes to SubstrATE will be documented here.

## [0.1.0] — unreleased

### Initial release

#### Core pipeline
- Full pipeline from genome FASTA to phylogenetic trees and synteny plots
- 26 built-in substrates covering marine, terrestrial, gut, and general
  polysaccharides
- Three PUL classification modes: bacteroidetes, generic, cazyme_only
- Survey mode for interactive substrate selection from dbCAN output
- Automatic gene prediction via Prodigal for nucleotide input
- dbCAN easy_substrate integration for CAZyme annotation and CGC
  prediction in one step

#### Sequence placement
- EPA-ng sequence placement onto reference trees of characterised CAZymes
- HMMER profile alignment (hmmbuild + hmmalign) for query-to-reference
  alignment
- Both .jplace (full placement data) and .newick (best placement)
  outputs per family
- Per-sequence placement confidence scores (likelihood weight ratios)
- Automatic fallback to de novo IQ-TREE2 tree building for families
  without reference trees or with too few reference sequences
- --skip_placement flag for de novo mode on all families

#### Reference sequence database
- build-reference-db subcommand: downloads characterised CAZyme sequences
  from CAZy via NCBI Entrez API with caching and resume support
- build-reference-trees subcommand: builds IQ-TREE2 reference trees with
  optional subsampling using proportional subfamily representation
- Global per-family reference trees shared across substrates

#### Activity annotation
- Activity patterns derived automatically from dbCAN
  fam-substrate-mapping.tsv
- Pattern review report written alongside output for manual curation
- Pattern overlap warnings between substrates with configurable threshold
- --overlap_threshold flag (default: 5, set to 0 to suppress)

#### Visualisation
- iTOL annotation files: sample colour strip, localisation colour strip,
  activity colour strip, branch colours, leaf symbols, display labels
- Two-phase colour configuration: assign once, edit and regenerate
- clinker synteny plots with top-level family colouring (GH16_3 -> GH16)
- Accessory enzyme class colours (sulfatase, deacetylase etc.)
- SusC/SusD transporter category colours

#### CLI
- substrate run — full pipeline with survey mode
- substrate annotate — dbCAN annotation only
- substrate classify — PUL classification only
- substrate tree — alignment, trimming and tree building only
- substrate visualise — iTOL annotation files only
- substrate synteny — clinker synteny plot only
- substrate survey — substrate survey from existing dbCAN output
- substrate list-substrates — list all built-in substrates
- substrate build-reference-db — download CAZy reference sequences
- substrate build-reference-trees — build IQ-TREE2 reference trees

### Known limitations
- Archaeal genomes not supported (CGC prediction requires Bacteroidetes-
  style PUL architecture)
- Activity patterns are auto-derived and should be reviewed against
  your dataset before publication
- --max_colours flag for datasets with >19 samples planned for v0.2.0
- Transporter phylogenetic tree planned for v0.2.0
- Dynamic substrate TSV replacing hardcoded FAMILY_MAP planned for v0.2.0
