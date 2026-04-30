# SubstrATE — Substrate Annotation Tool for Enzymes

SubstrATE is a command-line pipeline for identifying and annotating
polysaccharide utilisation loci (PULs) and CAZyme gene clusters (CGCs)
in bacterial genomes, with a focus on substrate specificity.

Given a set of genome assemblies or protein FASTAs and a target
polysaccharide substrate, SubstrATE:

1. Annotates CAZymes using [dbCAN](https://dbcan.readthedocs.io/)
2. Classifies CGCs as canonical PULs, non-canonical CGCs, or
   ungrouped CAZymes based on transporter gene co-localisation
3. Assigns enzymatic activity labels from EXPASY and dbCAN databases
4. Extracts substrate-relevant sequences and places them onto
   reference trees of characterised CAZymes using EPA-ng, or builds
   de novo trees using IQ-TREE2 where reference trees are unavailable
5. Generates [iTOL](https://itol.embl.de/) annotation files for tree
   visualisation
6. Produces [clinker](https://github.com/gamcil/clinker) synteny plots
   for PUL comparison

SubstrATE supports 26 built-in substrates covering marine algal
polysaccharides, plant cell wall components, gut glycans, and general
carbohydrates. Custom substrates can be derived automatically from the
dbCAN fam-substrate-mapping database at runtime.

---

## Table of contents

- [Installation](#installation)
- [Databases](#databases)
- [Reference sequence database](#reference-sequence-database)
- [Quick start](#quick-start)
- [Usage](#usage)
  - [Full pipeline](#full-pipeline)
  - [Survey mode](#survey-mode)
  - [Individual steps](#individual-steps)
- [Output structure](#output-structure)
- [Supported substrates](#supported-substrates)
- [PUL classification modes](#pul-classification-modes)
- [Activity patterns](#activity-patterns)
- [Comparison with SACCHARIS](#comparison-with-saccharis)
- [Citation](#citation)
- [License](#license)

---

## Installation

SubstrATE requires conda. If you do not have conda installed, download
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) first.

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/substrATE.git
cd substrATE

# Create the conda environment and install SubstrATE
conda env create -f environment.yml
conda activate substrATE

# Verify installation
substrate --help
```

### Dependencies

All dependencies are installed automatically by the conda environment.
The main external tools are:

| Tool | Version | Purpose |
|---|---|---|
| dbCAN | 5.2.8 | CAZyme annotation and CGC prediction |
| Prodigal | latest | Gene prediction (nucleotide input only) |
| MAFFT | 7.525 | Multiple sequence alignment |
| trimAl | 1.5 | Alignment trimming |
| IQ-TREE2 | 3.1.1 | Phylogenetic tree inference |
| EPA-ng | 0.3.8 | Evolutionary placement of sequences |
| HMMER | 3.4 | Profile alignment for EPA placement |
| clinker | 0.0.32 | Synteny plot generation |

Python dependencies: `biopython`, `pandas`, `click`, `requests`,
`beautifulsoup4`.

---

## Databases

SubstrATE requires the following database files:

### dbCAN database

Download and set up the dbCAN database following the
[dbCAN documentation](https://dbcan.readthedocs.io/en/latest/):

```bash
mkdir -p ~/db
cd ~/db
# Follow dbCAN database setup instructions
```

### EXPASY enzyme.dat

Download the EXPASY enzyme database:

```bash
wget https://ftp.expasy.org/databases/enzyme/enzyme.dat
```

### TCDB family definitions

Download the TCDB family definitions file:

```bash
wget https://www.tcdb.org/public/tcdb -O tc_family_definitions.tsv
```

---

## Reference sequence database

SubstrATE uses a database of characterised CAZyme sequences from CAZy
for EPA-ng sequence placement. This database is not included in the
repository and must be built once after installation.

### Step 1 — Download characterised sequences

```bash
substrate build-reference-db \
    --email your@email.com \
    --api_key YOUR_NCBI_API_KEY \
    --output substrate/data/reference_seqs
```

An NCBI API key is recommended for faster downloads (10 requests/second
vs 3). Register at https://www.ncbi.nlm.nih.gov/account/

This downloads characterised enzyme sequences from CAZy for all 26
built-in substrates and caches them locally. The download takes
approximately 30-60 minutes depending on your connection and NCBI load.

### Step 2 — Build reference trees

```bash
substrate build-reference-trees \
    --threads 8 \
    --max_seqs 200 \
    --output substrate/data/reference_trees
```

Builds IQ-TREE2 reference trees for each CAZyme family. Families with
fewer than 4 characterised sequences are skipped and will use de novo
tree building during pipeline runs.

The `--max_seqs` flag limits sequences per family tree using
proportional subfamily representation. This reduces build time without
losing coverage of important subfamilies. Recommended: 150-200 for
routine use, no limit for publication-quality trees.

Build time depends on family sizes — expect 4-12 hours for all
families with `--max_seqs 200`.

### Running without reference trees

If reference trees have not been built, use `--skip_placement` to fall
back to de novo IQ-TREE2 tree building for all families:

```bash
substrate run --skip_placement ...
```

---

## Quick start

### With existing dbCAN output

```bash
substrate run \
    --substrate laminarin \
    --dbcan_output /path/to/cgc_output/ \
    --db_dir ~/db \
    --expasy /path/to/enzyme.dat \
    --tcdb /path/to/tc_family_definitions.tsv \
    --output results/
```

### From genome assemblies

SubstrATE accepts nucleotide assemblies (`.fna`, `.fasta`, `.fa`) or
protein FASTAs (`.faa`). For nucleotide input, Prodigal is run
automatically:

```bash
substrate run \
    --substrate laminarin \
    --genomes /path/to/genomes/ \
    --db_dir ~/db \
    --expasy /path/to/enzyme.dat \
    --tcdb /path/to/tc_family_definitions.tsv \
    --output results/
```

### Multiple substrates

```bash
substrate run \
    --substrate laminarin \
    --substrate alginate \
    --substrate fucoidan \
    --dbcan_output /path/to/cgc_output/ \
    --db_dir ~/db \
    --expasy /path/to/enzyme.dat \
    --tcdb /path/to/tc_family_definitions.tsv \
    --output results/
```

---

## Usage

### Full pipeline

```
substrate run [OPTIONS]

Options:
  --substrate TEXT         Target substrate(s). Can be specified multiple
                           times. If omitted, runs survey mode.
  --genomes PATH           Path to genome FASTA directory
  --dbcan_output PATH      Path to existing dbCAN cgc_output directory
                           (skips annotation step)
  --db_dir PATH            Path to dbCAN database directory  [required]
  --expasy PATH            Path to EXPASY enzyme.dat  [required]
  --tcdb PATH              Path to TCDB tc_family_definitions.tsv  [required]
  --ref_metadata PATH      Path to reference sequence metadata TSV
  --ref_seqs PATH          Path to reference sequence FASTA directory
  --output PATH            Base output directory  [required]
  --threads INTEGER        Threads for MAFFT, IQ-TREE2 and EPA-ng
                           [default: 8]
  --pul_mode CHOICE        PUL classification mode  [default: bacteroidetes]
  --min_substrate_cazymes  Minimum substrate CAZymes per CGC  [default: 2]
  --skip_tree              Skip all tree building
  --skip_placement         Use de novo IQ-TREE2 for all families instead
                           of EPA placement (useful if reference trees
                           not yet built)
  --skip_clinker           Skip clinker synteny plot
  --force                  Overwrite existing output files
  --overlap_threshold INT  Pattern overlap warning threshold  [default: 5]
  --substrate_terms TEXT   Search terms for custom substrate derivation
```

### Survey mode

If `--substrate` is omitted, SubstrATE runs dbCAN annotation first
and then surveys all substrates in the dbCAN database, presenting
a ranked table of hits in your dataset:

```bash
substrate run \
    --genomes /path/to/genomes/ \
    --db_dir ~/db \
    --expasy /path/to/enzyme.dat \
    --tcdb /path/to/tc_family_definitions.tsv \
    --output results/
```

```
  Substrate              Genomes    Canonical PULs    Non-canonical
  ─────────────────────────────────────────────────────────────────
  laminarin              9          14                3
  alginate               6          9                 1
  fucoidan               3          4                 1
  ...

Which substrates would you like to analyse?
  all        — all substrates with hits
  min:N      — substrates with >= N canonical PULs
  a,b,c      — comma-separated substrate names

> laminarin, alginate
```

The survey can also be run independently on existing dbCAN output:

```bash
substrate survey \
    --dbcan_output /path/to/cgc_output/ \
    --db_dir ~/db \
    --output results/
```

### Individual steps

```bash
# Annotation only
substrate annotate --genomes /path/to/genomes/ --db_dir ~/db \
    --output results/

# Classification only
substrate classify --substrate laminarin \
    --dbcan_output /path/to/cgc_output/ --db_dir ~/db --output results/

# Tree building only (uses existing sequences/ directory)
substrate tree --substrate laminarin --output results/ --threads 8

# iTOL annotations only (regenerate after editing colour config)
substrate visualise --substrate laminarin --output results/

# Clinker synteny only
substrate synteny --substrate laminarin --output results/

# List built-in substrates
substrate list-substrates

# Survey existing dbCAN output
substrate survey --dbcan_output /path/to/cgc_output/ --db_dir ~/db
```

---

## Output structure

```
results/
├── logs/                          # Tool log files
├── prodigal/                      # Prodigal output (if nucleotide input)
├── cgc_output/                    # dbCAN output per genome
└── laminarin/                     # Per-substrate outputs
    ├── laminarin_family_hits.tsv      # CAZyme family hits with localisation
    ├── laminarin_substrate_hits.tsv   # Substrate prediction hits
    ├── laminarin_activity_annotated.tsv  # Activity annotations
    ├── laminarin_colour_config.tsv    # Editable colour assignments
    ├── laminarin_pattern_review.tsv   # Activity pattern review report
    ├── sequences/                     # Per-family FASTA files
    ├── alignments/                    # MAFFT alignments (de novo only)
    ├── trimmed/                       # trimAl trimmed alignments
    ├── trees/                         # IQ-TREE2 treefiles (de novo)
    ├── placements/                    # EPA-ng placement outputs
    │   └── GH16/
    │       ├── GH16.jplace            # Full placement data
    │       ├── GH16.newick            # Best placement as Newick
    │       └── GH16.placement_scores.tsv
    ├── itol_annotations/              # iTOL annotation files
    ├── genbank/                       # GenBank files per CGC
    └── clinker/                       # Clinker HTML and TSV outputs
```

### Key output files

**`{substrate}_family_hits.tsv`** — All CAZyme family hits for the
substrate, with columns including sample, gene ID, matched family,
CGC ID, localisation (canonical_PUL / non_canonical_CGC / outside_CGC),
and activity annotation.

**`{substrate}_activity_annotated.tsv`** — Activity-annotated hits
including reference sequences, used as input for iTOL and tree
interpretation.

**`{substrate}_colour_config.tsv`** — Colour assignments for samples,
activities, and localisations. Edit this file and rerun
`substrate visualise` to regenerate iTOL annotations with custom
colours without rebuilding the tree.

**`placements/{family}/{family}.jplace`** — Full EPA-ng placement data
in jplace format, compatible with iTOL's native placement visualisation.

**`placements/{family}/{family}.newick`** — Best placement positions
converted to Newick format for use with standard tree viewers.

**`placements/{family}/{family}.placement_scores.tsv`** — Per-sequence
placement confidence scores (likelihood weight ratios).

**`clinker/{substrate}_all_cgcs.html`** — Interactive clinker synteny
plot comparing all qualifying CGCs for the substrate.

---

## Supported substrates

SubstrATE includes 26 built-in substrates. Run `substrate list-substrates`
to see the full list with family counts.

| Category | Substrates |
|---|---|
| Marine/algal | laminarin, agar, carrageenan, alginate, fucoidan, ulvan, porphyran |
| Plant/terrestrial | xylan, arabinoxylan, pectin, chitin, cellulose, starch, beta_mannan, lichenan, mixed_linkage_glucan |
| Fructans | inulin, levan |
| Alpha-glucans | glycogen, pullulan |
| Host glycans | chondroitin_sulfate, heparan_sulfate, hyaluronic_acid |
| General | sucrose, xyloglucan, arabinogalactan |

### Custom substrates

For substrates not in the built-in list, SubstrATE derives CAZyme
families automatically from the dbCAN fam-substrate-mapping database:

```bash
substrate run \
    --substrate fucoidan \
    --substrate_terms "fucoidan,sulfated fucan" \
    ...
```

---

## PUL classification modes

SubstrATE supports three classification modes set with `--pul_mode`:

**`bacteroidetes`** (default) — requires SusC/SusD-type transporters
(TCDB families 1.B.14 and 8.A.46) co-located with substrate CAZymes.
Designed for Bacteroidetes PUL systems.

**`generic`** — requires any TC gene co-located with substrate CAZymes.
Suitable for non-Bacteroidetes bacteria with CAZyme gene clusters.

**`cazyme_only`** — classifies CGCs based purely on substrate CAZyme
count, with no transporter requirement. Useful for organisms where
transporter co-localisation is not expected.

> **Note:** SubstrATE is designed for bacteria with PUL-type CAZyme
> gene cluster systems. Archaeal genomes are not currently supported
> due to fundamental differences in CGC architecture. Use
> `--pul_mode cazyme_only` with caution for non-Bacteroidetes organisms
> and interpret results accordingly.

---

## Activity patterns

SubstrATE uses activity patterns to filter CGCs during GenBank file
generation, retaining only CGCs with a minimum number of genes with
substrate-relevant enzymatic activities.

Patterns are stored in `substrate/data/activity_patterns.tsv` and
are derived automatically from the dbCAN fam-substrate-mapping
database. All patterns are initially marked `reviewed=False`.

After running the pipeline on your dataset, inspect the
`{substrate}_pattern_review.tsv` report in each substrate output
directory. Once you are satisfied the patterns are appropriate for
your data, set `reviewed=True` in `activity_patterns.tsv` to suppress
the auto-derived patterns warning.

To regenerate patterns from a new version of the dbCAN database:

```bash
python scripts/generate_patterns.py \
    --fam_sub_map /path/to/fam-substrate-mapping.tsv \
    --output substrate/data/activity_patterns.tsv
```

### Pattern overlap warnings

SubstrATE warns when activity patterns for a substrate overlap
significantly with another substrate, since this may affect how
results should be interpreted:

```bash
# Suppress overlap warnings
substrate run --overlap_threshold 0 ...

# Increase sensitivity
substrate run --overlap_threshold 3 ...
```

---

## Comparison with SACCHARIS

[SACCHARIS](https://github.com/saccharis/SACCHARIS_2) is the closest
published tool to SubstrATE. The key differences are:

**Scope:** SACCHARIS works at the level of individual CAZyme families —
you provide sequences from one family and it builds an annotated tree.
SubstrATE works at the level of complete polysaccharide utilisation loci,
identifying which CGCs in a genome are involved in degrading a specific
substrate and comparing them across genomes.

**Input:** SACCHARIS requires sequences you provide manually. SubstrATE
starts from genome assemblies or protein FASTAs and handles the entire
workflow from dbCAN annotation through to visualisation.

**PUL context:** SACCHARIS has no concept of PULs or CGCs and does not
consider genomic organisation. SubstrATE's central output is the
classification of CGCs as canonical PULs, non-canonical CGCs, or
ungrouped CAZymes — the biologically meaningful unit for understanding
polysaccharide degradation in Bacteroidetes.

**Substrate focus:** SACCHARIS infers substrate specificity from tree
topology. SubstrATE takes the opposite approach — you specify the
substrate of interest and the tool identifies which genomic loci are
relevant to that substrate.

**Summary:** SACCHARIS answers *"What substrate does this CAZyme act
on?"* SubstrATE answers *"Which complete gene clusters in these genomes
are involved in degrading this substrate, and how do they compare?"*
The tools are complementary rather than competing.

---

## Citation

If you use SubstrATE in your research, please cite:

> [CITATION PLACEHOLDER — add when preprint/paper is available]

Please also cite the underlying tools:

- **dbCAN**: Zheng J et al. (2023) dbCAN3: automated CAZyme and
  substrate annotation. *Nucleic Acids Research*.
- **MAFFT**: Katoh K & Standley DM (2013) MAFFT multiple sequence
  alignment software. *Molecular Biology and Evolution*.
- **trimAl**: Capella-Gutierrez S et al. (2009) trimAl: a tool for
  automated alignment trimming. *Bioinformatics*.
- **IQ-TREE2**: Minh BQ et al. (2020) IQ-TREE 2: new models and
  methods for phylogenetic inference. *Molecular Biology and Evolution*.
- **EPA-ng**: Barbera P et al. (2019) EPA-ng: massively parallel
  evolutionary placement of genetic sequences. *Systematic Biology*.
- **HMMER**: Eddy SR (2011) Accelerated profile HMM searches.
  *PLOS Computational Biology*.
- **clinker**: Gilchrist CLM & Chooi YH (2021) clinker & clustermap.js.
  *Bioinformatics*.

---

## License

SubstrATE is released under the MIT License. See [LICENSE](LICENSE) for
details.

---

## Contributing

Bug reports, feature requests, and pull requests are welcome via the
[GitHub issue tracker](https://github.com/YOUR_USERNAME/substrATE/issues).

When adding a new built-in substrate, please:
1. Add entries to `FAMILY_MAP` and `SUBSTRATE_TERMS` in
   `substrate/classify_pul.py`
2. Add the substrate to `SUBSTRATE_SEARCH_CONFIG` in
   `scripts/generate_patterns.py`
3. Regenerate `activity_patterns.tsv` and verify patterns on a
   representative dataset before setting `reviewed=True`
4. Run `substrate build-reference-db` to download characterised
   sequences for the new substrate families
5. Rebuild reference trees with `substrate build-reference-trees`
