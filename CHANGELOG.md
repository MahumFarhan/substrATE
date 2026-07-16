# Changelog

All notable changes to SubstrATE will be documented here.

## [Unreleased]

### Planning
- v1.2.0 FAMILY_MAP plan updated: rather than dynamic derivation from
  fam-substrate-mapping.tsv at runtime (which would include promiscuous
  multi-substrate families and inflate canonical/non-canonical CGC calls),
  v1.2.0 will replace the hardcoded FAMILY_MAP with a curated
  family_map.tsv file (analogous to activity_patterns.tsv) with reviewed
  and excluded columns, loaded at runtime but editable without touching
  source code. Deliberate exclusions of promiscuous families are preserved
  as explicit excluded=True entries with reasons.

### Added
- `--ref_mode [diverse|relevant]` flag for `substrate run` — `diverse`
  (default) subsamples reference sequences for subfamily diversity;
  `relevant` keeps references phylogenetically closest to genomic
  sequences via MAFFT alignment + identity ranking
- `subsample_by_relevance()` in extract_seqs.py implementing relevant mode
- Reference substrate/activity labels in iTOL annotations now use an
  EC-number-derived controlled vocabulary (built from EXPASY enzyme.dat,
  231 ECs, 100% coverage) instead of the raw `substrate` metadata field,
  with keyword-based fallback for ambiguous ECs
- `load_ref_substrate_map()` in itol.py, plus Greek-symbol
  normalisation (alpha/beta/kappa/iota -> α/β/κ/ι) for reference labels
- Reference sequences now receive a distinct muted colour palette in
  iTOL activity strips, separate from genomic activity colours, to
  avoid palette collisions
- `reference_seqs_date.txt` written by `build-reference-db` on scrape;
  displayed by `check-db` as the reference sequence scrape date
- README: "How CGCs are detected" and "What canonical_PUL /
  non_canonical_CGC / outside_CGC mean" subsections under PUL
  classification modes
- README: "Reference sequence subsampling modes" subsection documenting
  `--ref_mode diverse|relevant`

### Fixed
- `clinker.py` and `genbank.py`: `prodoric` gene type (PRODORIC database
  transcriptional regulator cross-references) now labelled as `Regulator`
  with a dedicated colour (`#9970ab`) instead of falling through to
  `hypothetical` (grey). Affects all existing clinker synteny plots --
  regenerate outputs to see correct colouring.
- `clinker.py` and `genbank.py`: `Peptidase` gene type now labelled
  explicitly as `Peptidase` with a dedicated colour (`#5aae61`) instead
  of falling through to `hypothetical`. Same regeneration note applies.
- Raw GH/CBM family counts from `overview.tsv` must scan all four
  annotation columns (`dbCAN_hmm`, `dbCAN_sub`, `DIAMOND`,
  `Recommend Results`) to stay consistent with `family_hits.tsv`'s
  `all_annotations` logic. Scanning only `dbCAN_hmm`/`dbCAN_sub`
  undercounts genes whose family annotation appears only in `DIAMOND`
  (e.g. `CBM1+GH3` in `DIAMOND` with `-` in `dbCAN_hmm`/`dbCAN_sub`).
- `cli.py` `--ref_seqs`/`--ref_metadata` defaults now point to bundled
  reference data so reference sequences are appended automatically
  without requiring explicit flags — previously defaulted to `None`,
  causing reference sequences to be silently omitted from iTOL
  annotations despite appearing correctly in trees
- `run_dbcan.py`: dbCAN meta mode internally expects a `.fna` extension;
  inputs with other extensions (`.fasta`, `.fa`) are now symlinked to a
  `.fna` path before invocation, fixing `FileNotFoundError` on otherwise
  valid nucleotide FASTA input
- `FAMILY_MAP` (classify_pul.py) audited against dbCAN's authoritative
  `fam-substrate-mapping.tsv` and extended with substrate-specific
  families that were missing, fixing false `outside_CGC` classifications
  for genuine PULs:
    - glycogen/pullulan/starch: +30 alpha-glucan families (GH65, GH97,
      GH4, GH27, GH63, GH66, GH70, GH71, GH76, GH77, GH87, GH119,
      GH122, GH126, GH133, GH15, GH176, AA13, starch/glycogen-binding
      CBMs). Confirmed against a published PUL in *Muricauda* sp.
      MAR_2010_75 that was being missed (GH13+GH65 CGC, GH65 previously
      unmapped)
    - laminarin: +CBM102/103 (explicit laminarin-binding), CBM43/52/56,
      GH128/132/152/157/158/64/81 (beta-1,3-glucanase/laminarinase
      activity)
    - lichenan: +GH26/148/44/51/4 (mixed-linkage beta-1,3/1,4-glucanase
      activity)
    - arabinogalactan: +GH147, GH182, GH183
    - arabinoxylan: +CE15, GH121, GH146, GH93
    - beta_mannan: +GH164
    - cellulose: +CBM104
    - inulin: +CBM66, GH68
    - levan: +CBM38, CBM66, GH91
    - pectin: +CBM77, PL27, PL3
    - ulvan: +PL24, PL25, PL28, PL37
    - xylan: +CE16, GH159, GH6
    - xyloglucan: +GH26
    - alginate: +PL44

## [0.1.0] - 2025-05-21

### Added
- download-patterns command — downloads latest curated activity_patterns.tsv
  from GitHub releases, checks installed version, prompts for confirmation
  before overwriting; --force flag skips prompt
- activity_patterns_version.txt written to substrate/data/ on download,
  tracks installed release tag for version comparison
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
- add_fragments() in align.py using MAFFT --add
- place_sequences() in tree.py using IQ-TREE2 with reference tree as starting topology
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

## [0.0.1] - 2025-04-24

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

### Fixed (session 2)
- dbCAN now runs in meta mode for nucleotide input, enabling correct
  CGC prediction via pyrodigal (root cause of all CGCs showing as
  outside_CGC)
- MAFFT --addfragments replaced with --add (--addfragments caused
  segfaults on large reference alignments)
- IQ-TREE2 --tree-fix removed from place mode; reference tree now used
  as starting topology (-t), which does not require all sequences to be
  present in the tree
- Step counter 8/7 bug fixed (gbk_step was 6 instead of 5 in default
  run)
- Log directory now created before opening log file in trim.py and
  align.py
- S1 sulfatase colour correctly applied for bifunctional CAZyme+Sulfatase
  genes (e.g. GH16_3+Sulfatase|S1_7)
- TC family IDs now mapped to display categories in clinker.py
- SUBSTRATE_TERMS populated for all 25 substrates (previously only 3
  entries, causing KeyError in multi-substrate runs)

### Removed (session 2)
- check_prodigal() removed from run_dbcan.py (Prodigal no longer used)
- run_prodigal() removed from run_dbcan.py (dead code)

### Tests (session 2)
- test_activity.py written from scratch: normalise_activity,
  extract_primary_ec, get_activity_label, load_ec_names,
  load_family_activities, annotate_references, load_patterns
- test_classify_pul.py: added TestSubstrateTerms and TestFamilyMap
  covering all 25 substrates and removal of mixed_linkage_glucan
- test_extract_seqs.py: no changes needed (signatures unchanged)

### Fixed (session 3)
- check_prodigal and needs_prodigal references removed from cli.py
  _validate_tools() — caused AttributeError for any nucleotide input
  after check_prodigal was removed from run_dbcan.py in session 2

### Curation (session 3)
- Batch 1 activity patterns curated: starch, cellulose, chitin, xylan,
  pectin, beta_mannan, arabinoxylan, xyloglucan
- 84 patterns marked reviewed=True across 8 substrates
- Strict/permissive mode assigned per pattern based on substrate
  specificity
- Reviewed pattern files archived in curation/activity_patterns/
