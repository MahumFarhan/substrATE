# Changelog

All notable changes to SubstrATE will be documented here.

## [Unreleased]

### Added
- Strict/permissive activity pattern modes via --pattern_mode flag
  (default: permissive)
- load_patterns() as single authoritative pattern loader in
  parse_substrates.py, replacing direct CSV reads across modules
- Auto-detection of tree building mode per CAZyme family:
    - place: MAFFT --addfragments + IQ-TREE2 --tree-fix onto reference
      tree (fastest, reproducible across runs)
    - merge: de novo with CAZy reference sequences merged in (fallback
      when no reference tree exists)
    - denovo: genomic sequences only (fallback or via --denovo flag)
- add_fragments() in align.py using MAFFT --addfragments
- place_sequences() in tree.py using IQ-TREE2 --tree-fix
- --denovo flag to force de novo tree building for all families
- normalise_db_dir() in run_dbcan.py to handle case mismatches between
  dbCAN expected filenames and Windows/WSL2 filesystem conventions
- mode column in activity_patterns.tsv (strict/permissive per token)
- pattern_mode column in pattern review reports
- annotate_references() now supports both legacy and current
  reference_metadata.tsv schemas (protein_name + ec_numbers)

### Changed
- Tree building now auto-selects place/merge/denovo per family with
  clear per-family messages showing mode and sequence counts
- --tree_mode flag removed (replaced by auto-detection + --denovo)
- SusC and TonB-dep. colour entries merged in default_colours.tsv
- Pattern overlap warnings now report the active pattern_mode

### Removed
- mixed_linkage_glucan removed from FAMILY_MAP and activity_patterns.tsv
  (patterns too generic for reliable classification)
- place_sequences.py removed (EPA-ng based, superseded by tree.py
  IQ-TREE2 --tree-fix implementation)
- epa-ng and hmmer removed from environment.yml dependencies

### Fixed
- annotate_references() KeyError when activity column absent in new
  reference_metadata.tsv schema
- dbCAN database filename case mismatches on WSL2/Windows filesystems
  (TCDB.dmnd, CAZy.dmnd, dbCAN.hmm, dbCAN-sub.hmm)
- Noise tokens removed from activity_patterns.tsv across all substrates
- Parse artifacts removed (glycosaminoglycan, with trailing comma)

## [0.1.0] - 2025-04-24

### Initial release
- Full pipeline from genome FASTA to phylogenetic trees and synteny plots
- 25 built-in substrates covering marine, terrestrial, gut, and general
  polysaccharides
- Three PUL classification modes: bacteroidetes, generic, cazyme_only
- Survey mode for interactive substrate selection from dbCAN output
- Automatic gene prediction via Prodigal for nucleotide input
- dbCAN easy_substrate integration for CAZyme annotation and CGC
  prediction in one step
- De novo IQ-TREE2 trees per CAZyme family
- build-reference-db subcommand: downloads characterised CAZyme sequences
  from CAZy via NCBI Entrez API
- build-reference-trees subcommand: builds IQ-TREE2 reference trees with
  proportional subfamily subsampling
- Activity patterns auto-derived from dbCAN fam-substrate-mapping.tsv
- Pattern review reports written alongside output for manual curation
- Pattern overlap warnings with configurable threshold
- iTOL annotation files: sample colours, localisation, activity,
  branch colours, leaf symbols, display labels
- clinker synteny plots with top-level family colouring
- Accessory enzyme class colours (sulfatase, deacetylase etc.)
- SusC/SusD transporter category colours
- substrate run, annotate, classify, tree, visualise, synteny,
  survey, list-substrates, build-reference-db, build-reference-trees

### Known limitations
- Archaeal genomes not supported
- Activity patterns should be reviewed against your dataset before
  publication
- --max_colours flag for datasets with >19 samples planned for v0.2.0
- Dynamic substrate TSV replacing hardcoded FAMILY_MAP planned for v0.2.0
- Transporter phylogenetic tree planned for v0.2.0
