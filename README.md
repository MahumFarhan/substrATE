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
4. Extracts substrate-relevant sequences, merges them with
   characterised CAZyme reference sequences from CAZy, and builds
   phylogenetic trees using IQ-TREE2
5. Generates [iTOL](https://itol.embl.de/) annotation files for tree
   visualisation
6. Produces [clinker](https://github.com/gamcil/clinker) synteny plots
   for PUL comparison

SubstrATE supports 25 built-in substrates covering marine algal
polysaccharides, plant cell wall components, gut glycans, and general
carbohydrates. Custom substrates can be derived automatically from the
dbCAN fam-substrate-mapping database at runtime.

> **Maintenance notice:** SubstrATE is developed and maintained on a
> best-effort basis. Issues and pull requests may not always receive a
> timely response. The tool is provided as-is for the community to use
> and adapt freely under the MIT licence.

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
- [Acknowledgements](#acknowledgements)
- [License](#license)

---

## Installation

SubstrATE requires conda. If you do not have conda installed, download
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) first.

```bash
# Clone the repository
git clone https://github.com/MahumFarhan/substrATE.git
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
| MAFFT | 7.525 | Multiple sequence alignment |
| trimAl | 1.5 | Alignment trimming |
| IQ-TREE2 | 3.1.1 | Phylogenetic tree inference |
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
to enrich phylogenetic trees with biological context. This database is
not included in the repository and must be built once after
installation.

### Step 1 — Download characterised sequences

```bash
substrate build-reference-db \
    --email your@email.com \
    --api_key YOUR_NCBI_API_KEY \
    --output substrate/data/reference_seqs
```

An NCBI API key is recommended for faster downloads (10 requests/second
vs 3). Register at https://www.ncbi.nlm.nih.gov/account/

This downloads characterised enzyme sequences from CAZy for all 25
built-in substrates. The download takes approximately 30–60 minutes
depending on your connection and NCBI load.

### Step 2 — Build reference trees (optional)

Reference trees are pre-built from characterised sequences and used
when SubstrATE auto-selects place mode for a CAZyme family. Place mode
fixes the reference topology and places genomic sequences onto it,
giving reproducible trees that are directly comparable across runs.

```bash
substrate build-reference-trees \
    --threads 8 \
    --max_seqs 150 \
    --output substrate/data/reference_trees
```

This step is not required if you are happy with the default merge
behaviour, where SubstrATE builds a de novo tree per family with
reference sequences merged in automatically.

### Tree building modes

SubstrATE auto-selects a tree building mode per CAZyme family based
on what reference data is available:

| Mode | Description | When selected |
|---|---|---|
| `place` | Genomic sequences placed onto fixed reference tree topology | Reference tree exists for family |
| `merge` | De novo tree with reference sequences merged in | Reference seqs exist, no reference tree |
| `denovo` | De novo tree from genomic sequences only | No reference data, or --denovo flag set |

Use `--denovo` to force de novo building for all families regardless
of available reference data.

> **Note on tree quality:** Reference trees are built with
> `LG+G4 --fast` (no bootstrap) — sufficient as a placement backbone
> but not publication quality. Per-family pipeline trees use full model
> selection and bootstrap (`-m TEST -B 1000`).
>
> **Note on reproducibility:** Merge and denovo mode trees are not
> fully reproducible between runs — IQ-TREE2 uses stochastic
> hill-climbing and random starting trees. For reproducible results,
> use `--seed` to fix the random seed:
> ```bash
> substrate run --seed 42 ...
> ```
> Place mode trees are fully reproducible without a seed since the
> reference topology is fixed.

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
protein FASTAs (`.faa`). For nucleotide input, dbCAN runs in meta mode
and handles gene prediction internally using pyrodigal:

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
--threads INTEGER        Threads for MAFFT and IQ-TREE2  [default: 8]
--pul_mode CHOICE        PUL classification mode  [default: bacteroidetes]
--min_substrate_cazymes  Minimum substrate CAZymes per CGC  [default: 2]
--pattern_mode CHOICE    Activity pattern stringency  [default: permissive]
--skip_tree              Skip alignment, trimming, and tree building
--skip_clinker           Skip clinker synteny plot
--denovo                 Force de novo tree building for all families
--force                  Overwrite existing output files
--overlap_threshold INT  Pattern overlap warning threshold  [default: 5]
--substrate_terms TEXT   Search terms for custom substrate derivation
--max_colours INT        Maximum samples before HSL colour generation

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

laminarin, alginate


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

# List built-in substrates and reference sequence counts
substrate list-substrates
substrate family-sizes --substrate laminarin

# Test installation on bundled Gramella forsetii genome
substrate test-install

# Survey existing dbCAN output
substrate survey --dbcan_output /path/to/cgc_output/ --db_dir ~/db
```

---

## Output structure
results/
├── logs/                          # Tool log files
├── cgc_output/                    # dbCAN output per genome
└── laminarin/                     # Per-substrate outputs
├── laminarin_family_hits.tsv      # CAZyme family hits with localisation
├── laminarin_substrate_hits.tsv   # Substrate prediction hits
├── laminarin_activity_annotated.tsv  # Activity annotations
├── laminarin_colour_config.tsv    # Editable colour assignments
├── laminarin_pattern_review.tsv   # Activity pattern review report
├── sequences/                     # Per-family FASTA files
├── alignments/                    # MAFFT alignments
├── trimmed/                       # trimAl trimmed alignments
├── trees/                         # IQ-TREE2 treefiles
├── itol_annotations/              # iTOL annotation files
├── genbank/                       # GenBank files per CGC
└── clinker/                       # Clinker HTML and TSV outputs

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
colours without rebuilding trees.

**`clinker/{substrate}_all_cgcs.html`** — Interactive clinker synteny
plot comparing all qualifying CGCs for the substrate.

---

## Supported substrates

SubstrATE includes 25 built-in substrates. Run `substrate list-substrates`
to see the full list with family counts.

| Category | Substrates |
|---|---|
| Marine/algal | laminarin, agar, carrageenan, alginate, fucoidan, ulvan, porphyran |
| Plant/terrestrial | xylan, arabinoxylan, pectin, chitin, cellulose, starch, beta_mannan, lichenan, xyloglucan |
| Fructans | inulin, levan |
| Alpha-glucans | glycogen, pullulan |
| Host glycans | chondroitin_sulfate, heparan_sulfate, hyaluronic_acid |
| General | sucrose, arabinogalactan |

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
> gene cluster systems. Archaeal genomes are not currently supported.
> Use `--pul_mode cazyme_only` with caution for non-Bacteroidetes
> organisms and interpret results accordingly.
>
> ### Metagenomic input

SubstrATE supports metagenomic-assembled genomes (MAGs) as input, with
the following expectations and caveats:

- **Binned MAGs only.** Each input file must be a single-organism bin.
  Unbinned metagenomic assemblies are not supported — CGCFinder links
  genes into clusters based on genomic proximity, and fragmented or
  mixed-organism contigs will produce incomplete or spurious CGCs.

- **Use `--pul_mode generic` or `--pul_mode cazyme_only`.** The default
  `bacteroidetes` mode requires SusC/SusD-type transporters and will
  miss PULs from non-Bacteroidetes organisms common in metagenomic
  datasets.

- **Contig fragmentation reduces CGC recovery.** PULs that span contig
  boundaries will be split or missed entirely. This is a fundamental
  limitation of short-read assemblies and cannot be resolved within
  substrATE.

- **Taxonomy is not built in.** To link substrATE results to taxonomic
  annotations, join your taxonomy table against
  `{substrate}_family_hits.tsv` on the sample name column. Tools such
  as GTDB-Tk are recommended for MAG taxonomic classification.

---

## Activity patterns

SubstrATE uses activity patterns to filter CGCs during GenBank file
generation, retaining only CGCs with a minimum number of genes bearing
substrate-relevant enzymatic activities.

Patterns are stored in `substrate/data/activity_patterns.tsv` and
are derived automatically from the dbCAN fam-substrate-mapping
database. All patterns are initially marked `reviewed=False`.

After running the pipeline, inspect the
`{substrate}_pattern_review.tsv` report in each substrate output
directory. Once satisfied, set `reviewed=True` in
`activity_patterns.tsv` to suppress the auto-derived patterns warning.

### Pattern stringency

Activity patterns are assigned one of two modes, which determine how
they are applied during CGC filtering:

**`permissive`** — the pattern is relevant for the substrate but the
enzyme may also act on other substrates. Applied in both permissive
and strict pipeline runs. Examples: `arabinosidase` for arabinoxylan
(also acts on arabinogalactan), `glucosidase` for starch (broad
substrate range).

**`strict`** — the pattern is highly specific to the substrate and
unlikely to produce false positives. Only applied when running with
`--pattern_mode strict`. Examples: `xylanase` for xylan,
`laminarinase` for laminarin, `chitinase` for chitin.

Use `--pattern_mode strict` for conservative filtering that reduces
false positives from enzymes active on multiple substrates. Use the
default `--pattern_mode permissive` for maximum sensitivity.

```bash
substrate run --pattern_mode strict ...
```

To regenerate patterns from a new version of the dbCAN database:

```bash
python scripts/generate_patterns.py \
    --fam_sub_map /path/to/fam-substrate-mapping.tsv \
    --output substrate/data/activity_patterns.tsv
```

### Pattern overlap warnings

SubstrATE warns when activity patterns for a substrate overlap
significantly with another substrate:

```bash
substrate run --overlap_threshold 0 ...   # suppress warnings
substrate run --overlap_threshold 3 ...   # increase sensitivity
```

---

## Acknowledgements

Please cite the underlying tools that SubstrATE depends on:

- **dbCAN**: Zheng J et al. (2023) dbCAN3: automated CAZyme and
  substrate annotation. *Nucleic Acids Research*.
- **MAFFT**: Katoh K & Standley DM (2013) MAFFT multiple sequence
  alignment software. *Molecular Biology and Evolution*.
- **trimAl**: Capella-Gutierrez S et al. (2009) trimAl: a tool for
  automated alignment trimming. *Bioinformatics*.
- **IQ-TREE2**: Minh BQ et al. (2020) IQ-TREE 2: new models and
  methods for phylogenetic inference. *Molecular Biology and Evolution*.
- **clinker**: Gilchrist CLM & Chooi YH (2021) clinker & clustermap.js.
  *Bioinformatics*.

---

## License

SubstrATE is released under the MIT License. See [LICENSE](LICENSE) for
details.

---

## Contributing

Bug reports and pull requests are welcome via the
[GitHub issue tracker](https://github.com/MahumFarhan/substrATE/issues).
Please note that responses may be delayed and fixes are not guaranteed.

When adding a new built-in substrate:
1. Add entries to `FAMILY_MAP` and `SUBSTRATE_TERMS` in
   `substrate/classify_pul.py`
2. Add the substrate to `SUBSTRATE_SEARCH_CONFIG` in
   `scripts/generate_patterns.py`
3. Regenerate `activity_patterns.tsv` and verify patterns on a
   representative dataset before setting `reviewed=True`
4. Run `substrate build-reference-db` to download characterised
   sequences for the new substrate families
5. Rebuild reference trees with `substrate build-reference-trees`
