# Activity Pattern Curation

This directory documents the curation of activity patterns used by
SubstrATE to filter CAZyme gene clusters (CGCs) during analysis.

Activity patterns are short strings matched against enzymatic activity
labels assigned to each gene. They are used to retain only CGCs with a
minimum number of genes bearing substrate-relevant enzymatic activities,
reducing false positives from non-specific CAZyme family hits.

---

## What is an activity pattern?

Each pattern is a case-insensitive substring matched against the
`activity` column in `{substrate}_activity_annotated.tsv`. For example,
the pattern `laminarinase` matches any gene annotated as
`laminarinase`, `endo-laminarinase`, `exo-laminarinase`, etc.

Patterns are stored in `substrate/data/activity_patterns.tsv` with the
following columns:

| Column | Description |
|---|---|
| `substrate` | Substrate name (must match a key in `FAMILY_MAP`) |
| `pattern` | Substring to match against activity labels |
| `source` | `auto_derived` or `curated` |
| `reviewed` | `True` if manually reviewed, `False` if auto-derived |
| `mode` | `permissive` (default) or `strict` |

Patterns with `mode=strict` are only applied when the pipeline is run
with `--pattern_mode strict`, which reduces false positives at the cost
of sensitivity.

---

## Curation workflow

### Step 1 — Download reference strain genomes

Reference strains with well-characterised polysaccharide degradation
systems are used to generate and validate patterns. The full list of
recommended strains with NCBI accessions is in `reference_strains.txt`.

Download genomes using the NCBI datasets tool:

```bash
./datasets download genome accession GCA_XXXXXXXXX.X --include genome
```

Protein FASTAs (`.faa`) can be used directly if available, skipping
gene prediction.

### Step 2 — Run the pipeline on reference strains

Run SubstrATE on the reference genome(s) for the substrate being
curated:

```bash
substrate run \
    --substrate xylan \
    --genomes /path/to/reference/genomes/ \
    --db_dir /path/to/dbcan/db \
    --expasy /path/to/enzyme.dat \
    --tcdb /path/to/tc_family_definitions.tsv \
    --output /path/to/curation_output/
```

### Step 3 — Review the pattern report

Each run produces a `{substrate}_pattern_review.tsv` in the substrate
output directory. This report lists every pattern, how many genes it
matched, and the activity labels of those genes. Open it in a
spreadsheet application and review each pattern:

- Does the matched activity make sense for this substrate?
- Is the pattern too broad (matching unrelated activities)?
- Is the pattern too narrow (missing known activities)?

Cross-check against the primary literature for the reference strain
where possible.

### Step 4 — Update activity_patterns.tsv

For patterns that are correct and specific:
- Set `reviewed=True`
- Set `mode=strict` if the pattern is highly specific and unlikely to
  produce false positives across other substrates
- Leave `mode=permissive` for broader patterns that are useful but
  may overlap with other substrates

For patterns that are incorrect or too broad:
- Remove them from `activity_patterns.tsv`

For activities present in the reference strain output but missing from
the pattern list:
- Add them manually with `source=curated` and `reviewed=True`

### Step 5 — Commit after each batch

Commit `activity_patterns.tsv` after each substrate batch so progress
is tracked and reversible:

```bash
git add substrate/data/activity_patterns.tsv
git commit -m "curation: batch 1 — xylan, arabinoxylan, pectin, starch"
```

---

## Curation batches

Curation is organised in three batches, prioritised by substrate
importance and reference strain availability.

### Batch 1 — Plant/terrestrial substrates
Substrates: starch, cellulose, chitin, xylan, pectin, beta_mannan,
arabinoxylan, xyloglucan

Reference strains used:
- *Bacteroides thetaiotaomicron* VPI-5482 (GCA_000011065.1) — starch, pectin
- *Bacteroides ovatus* ATCC 8483 (GCA_000310825.1) — xylan, arabinoxylan, beta_mannan, cellulose, starch
- *Flavobacterium johnsoniae* UW101 (GCA_000016645.1) — chitin
- *Cytophaga hutchinsonii* ATCC 33406 (GCA_000015065.1) — cellulose
- *Bacteroides cellulosilyticus* WH2 (GCA_000373585.1) — xyloglucan

### Batch 2 — Marine/algal and gut substrates
Substrates: laminarin, alginate, carrageenan, fucoidan, heparan_sulfate,
chondroitin_sulfate, hyaluronic_acid, inulin, levan, pullulan

Reference strains used:
- *Zobellia galactanivorans* DsiJT (GCA_000227045.1) — laminarin, alginate, carrageenan, fucoidan
- *Bacteroides thetaiotaomicron* VPI-5482 (GCA_000011065.1) — heparan_sulfate, chondroitin_sulfate, hyaluronic_acid, inulin, levan
- *Wenyingzhuangia fucanilytica* CZ1127 (GCA_001682895.1) — fucoidan
- *Formosa agariphila* KMM 3901 (GCA_000512255.1) — ulvan (Batch 3 overlap)

### Batch 3 — Remaining substrates
Substrates: agar, porphyran, ulvan, lichenan, arabinogalactan,
glycogen, sucrose

Reference strains used:
- *Zobellia galactanivorans* DsiJT (GCA_000227045.1) — agar, porphyran
- *Formosa agariphila* KMM 3901 (GCA_000512255.1) — ulvan
- *Bacteroides ovatus* ATCC 8483 (GCA_000310825.1) — lichenan
- *Bacteroides thetaiotaomicron* VPI-5482 (GCA_000011065.1) — glycogen, sucrose

---

## Downloading curated patterns

A curated `activity_patterns.tsv` covering all 25 built-in substrates
is available for download without running curation yourself:

```bash
substrate download-patterns
```

This command downloads the reviewed pattern file from the latest
SubstrATE GitHub release and places it in `substrate/data/`. See the
main README for details.

> **Note:** `substrate download-patterns` is not yet implemented.
> It will be added in a future release once curation is complete.
> In the meantime, the curated file can be downloaded manually from
> the [GitHub releases page](https://github.com/MahumFarhan/substrATE/releases).
