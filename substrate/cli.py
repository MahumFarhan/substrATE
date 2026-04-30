"""
substrATE command line interface.

Entry point for all pipeline subcommands. The main `run` command
executes the full pipeline for one or more substrates. Individual
subcommands allow running specific steps independently using
intermediate files from a previous run.
"""
import os
import sys
import traceback
import pandas as pd
import click

from substrate import (run_dbcan, classify_pul, activity,
                       extract_seqs, align, trim, tree,
                       genbank, itol, clinker, parse_substrates,
                       place_sequences)
from substrate.align import TooFewSequencesError
from substrate.align import ToolNotFoundError as AlignToolError
from substrate.trim  import ToolNotFoundError as TrimToolError
from substrate.tree  import ToolNotFoundError as TreeToolError
from substrate.clinker import ToolNotFoundError as ClinkerToolError
from substrate.run_dbcan import ToolNotFoundError as DbcanToolError
from substrate.place_sequences import (ToolNotFoundError as EpaToolError,
                                        place_all_families)


# ── Constants ─────────────────────────────────────────────────────────────────

_DATA_DIR     = os.path.join(os.path.dirname(__file__), 'data')
_COLOURS_FILE = os.path.join(_DATA_DIR, 'default_colours.tsv')
_PATTERNS_FILE = os.path.join(_DATA_DIR, 'activity_patterns.tsv')


# ── Helpers ───────────────────────────────────────────────────────────────────

def _section(title, width=60):
    """Print a section header."""
    click.echo(f"\n{'='*width}")
    click.echo(f"  {title}")
    click.echo(f"{'='*width}")


def _step(n, total, description):
    """Print a step header."""
    click.echo(f"\nStep {n}/{total}: {description}...")


def _success(msg):
    click.echo(f"  ✓ {msg}")


def _warn(msg):
    click.echo(f"  WARNING: {msg}", err=True)


def _error(msg):
    click.echo(f"  ERROR: {msg}", err=True)


def _validate_paths(**kwargs):
    """
    Validate that required paths exist.

    Args:
        **kwargs: name=path pairs to check

    Raises:
        click.ClickException if any path does not exist
    """
    for name, path in kwargs.items():
        if path and not os.path.exists(path):
            raise click.ClickException(
                f"--{name} path does not exist: {path}"
            )


def _validate_tools(skip_tree=False, skip_clinker=False,
                    needs_prodigal=False):
    """
    Check all required external tools are available.

    Args:
        skip_tree:      if True, skip MAFFT/trimAl/IQ-TREE2 checks
        skip_clinker:   if True, skip clinker check
        needs_prodigal: if True, check for Prodigal

    Raises:
        click.ClickException if any required tool is missing
    """
    click.echo("\nChecking required tools...")

    checks = [
        ('dbCAN',   run_dbcan.check_dbcan),
    ]
    if needs_prodigal:
        checks.append(('Prodigal', run_dbcan.check_prodigal))
    if not skip_tree:
        checks.extend([
            ('MAFFT',   align.check_mafft),
            ('trimAl',  trim.check_trimal),
            ('IQ-TREE2', lambda: tree.check_iqtree()[1]),
        ])
    if not skip_clinker:
        checks.append(('clinker', clinker.check_clinker))

    all_ok = True
    for name, check_fn in checks:
        try:
            version = check_fn()
            click.echo(f"  ✓ {name}: {version}")
        except Exception as e:
            _error(f"{name} not found: {e}")
            all_ok = False

    if not all_ok:
        raise click.ClickException(
            "One or more required tools are missing. "
            "Install them with conda env create -f environment.yml"
        )


def _parse_substrate_selection(response, available):
    """
    Parse user substrate selection from survey prompt.

    Accepts:
      'all'      — all substrates with hits
      'min:N'    — all substrates with >= N canonical PULs
      'a,b,c'    — comma-separated substrate names

    Args:
        response:  user input string
        available: DataFrame from survey_substrates()

    Returns:
        list of selected substrate name strings

    Raises:
        click.ClickException if input is invalid
    """
    response = response.strip().lower()

    if response == 'all':
        return available['substrate'].tolist()

    if response.startswith('min:'):
        try:
            threshold = int(response.split(':')[1])
            selected = available[
                available['n_canonical_pul'] >= threshold
            ]['substrate'].tolist()
            if not selected:
                raise click.ClickException(
                    f"No substrates have >= {threshold} canonical PULs"
                )
            return selected
        except (ValueError, IndexError):
            raise click.ClickException(
                f"Invalid min: format. Use e.g. 'min:3'"
            )

    # Comma-separated names
    selected = [s.strip() for s in response.split(',')
                if s.strip()]
    available_names = set(available['substrate'].tolist())
    invalid = [s for s in selected if s not in available_names]
    if invalid:
        raise click.ClickException(
            f"Unknown substrate(s): {', '.join(invalid)}\n"
            f"Available: {', '.join(sorted(available_names))}"
        )
    return selected


def _get_fam_sub_map(db_dir):
    """Find fam-substrate-mapping.tsv in db_dir."""
    candidates = [
        os.path.join(db_dir, 'fam-substrate-mapping.tsv'),
        os.path.join(db_dir, 'fam_substrate_mapping.tsv'),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    raise click.ClickException(
        f"Cannot find fam-substrate-mapping.tsv in {db_dir}. "
        f"Check your --db_dir path."
    )


def _save_dataframe(df, path, description):
    """Save a DataFrame to TSV with a status message."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, sep='\t', index=False)
    _success(f"Saved {description}: {path}")


# ── run command ───────────────────────────────────────────────────────────────

@click.group()
def main():
    """substrATE: Substrate Annotation Tool for Enzymes."""
    pass


@main.command()
@click.option('--substrate', multiple=True,
              help='Target substrate(s). Can be specified multiple times. '
                   'If omitted, dbCAN annotation runs first and a survey '
                   'prompts you to select substrates interactively.')
@click.option('--genomes', type=click.Path(), default=None,
              help='Path to directory containing genome FASTA files')
@click.option('--dbcan_output', type=click.Path(), default=None,
              help='Path to existing dbCAN cgc_output directory '
                   '(skips annotation step)')
@click.option('--db_dir', required=True, type=click.Path(),
              help='Path to dbCAN database directory')
@click.option('--expasy', required=True, type=click.Path(),
              help='Path to EXPASY enzyme.dat file')
@click.option('--tcdb', required=True, type=click.Path(),
              help='Path to TCDB tc_family_definitions.tsv file')
@click.option('--ref_metadata', type=click.Path(), default=None,
              help='Path to reference sequence metadata TSV')
@click.option('--ref_seqs', type=click.Path(), default=None,
              help='Path to directory containing characterised '
                   'reference sequences')
@click.option('--output', required=True, type=click.Path(),
              help='Base output directory. Each substrate gets its own '
                   'subdirectory: <output>/<substrate>/')
@click.option('--threads', default=8, show_default=True,
              help='Number of threads for MAFFT and IQ-TREE2')
@click.option('--pul_mode', default='bacteroidetes', show_default=True,
              type=click.Choice(['bacteroidetes', 'generic',
                                 'cazyme_only']),
              help='PUL classification mode')
@click.option('--min_substrate_cazymes', default=2, show_default=True,
              help='Minimum substrate CAZymes required per CGC')
@click.option('--skip_tree', is_flag=True, default=False,
              help='Skip all tree building (placement and de novo)')
@click.option('--skip_placement', is_flag=True, default=False,
              help='Skip EPA placement, use de novo IQ-TREE2 for all '
                   'families (useful if reference trees not yet built)')
@click.option('--skip_clinker', is_flag=True, default=False,
              help='Skip clinker synteny plot')
@click.option('--force', is_flag=True, default=False,
              help='Overwrite existing output files')
@click.option('--overlap_threshold', default=5, show_default=True,
              help='Minimum shared activity patterns to trigger overlap '
                   'warning. Set to 0 to suppress.')
@click.option('--substrate_terms', default=None,
              help='Comma-separated search terms for custom substrate '
                   'family derivation (only needed for substrates not '
                   'in the built-in list)')
def run(substrate, genomes, dbcan_output, db_dir, expasy, tcdb,
        ref_metadata, ref_seqs, output, threads, pul_mode,
        min_substrate_cazymes, skip_tree, skip_placement,
        skip_clinker, force, overlap_threshold, substrate_terms):
    """Run the full analysis pipeline."""

    os.makedirs(output, exist_ok=True)

    # ── Validate input paths ──────────────────────────────────────────────────
    _validate_paths(
        db_dir=db_dir,
        expasy=expasy,
        tcdb=tcdb,
        ref_metadata=ref_metadata,
        ref_seqs=ref_seqs,
        dbcan_output=dbcan_output,
        genomes=genomes,
    )

    # Locate fam-substrate-mapping.tsv
    fam_sub_map = _get_fam_sub_map(db_dir)

    # Determine if Prodigal is needed
    needs_prodigal = (
        genomes is not None and
        dbcan_output is None and
        any(f.endswith(('.fna', '.fasta', '.fa'))
            for f in os.listdir(genomes)
            if os.path.isfile(os.path.join(genomes, f)))
    )

    # ── Validate tools ────────────────────────────────────────────────────────
    _validate_tools(
        skip_tree=skip_tree,
        skip_clinker=skip_clinker,
        needs_prodigal=needs_prodigal,
    )

    # ── dbCAN annotation ──────────────────────────────────────────────────────
    if dbcan_output:
        cgc_output_dir = dbcan_output
        _success(f"Using existing dbCAN output: {cgc_output_dir}")
    elif genomes:
        _section("Step 1: dbCAN Annotation")
        cgc_output_dir = run_dbcan.annotate_genomes(
            input_dir=genomes,
            output_dir=output,
            db_dir=db_dir,
            threads=threads,
            force=force,
        )
    else:
        raise click.ClickException(
            "Either --genomes or --dbcan_output must be provided."
        )

    # ── Survey and substrate selection ────────────────────────────────────────
    if not substrate:
        _section("Step 2: Substrate Survey")
        click.echo("No substrates specified — surveying dbCAN output...\n")

        survey_df = parse_substrates.survey_substrates(
            cgc_output_dir=cgc_output_dir,
            fam_sub_map=fam_sub_map,
            pul_mode=pul_mode,
            min_cazymes=min_substrate_cazymes,
        )

        total = len(pd.read_csv(
            fam_sub_map, sep='\t',
            usecols=[1]).iloc[:, 0].dropna().unique())

        parse_substrates.print_survey_results(survey_df, total)

        if survey_df.empty:
            raise click.ClickException(
                "No substrate hits found in dbCAN output. "
                "Check your input files and --db_dir path."
            )

        # Save survey results
        survey_path = os.path.join(output, 'survey_results.tsv')
        survey_df.to_csv(survey_path, sep='\t', index=False)
        _success(f"Survey results saved to {survey_path}")

        # Interactive substrate selection
        click.echo(
            "\nWhich substrates would you like to analyse?\n"
            "Options:\n"
            "  all        — all substrates with hits\n"
            "  min:N      — substrates with >= N canonical PULs\n"
            "  a,b,c      — comma-separated substrate names\n"
        )
        response = click.prompt("Your selection")
        substrate = _parse_substrate_selection(response, survey_df)
        click.echo(f"\nSelected: {', '.join(substrate)}")

    # ── Pattern overlap check ─────────────────────────────────────────────────
    if overlap_threshold > 0:
        parse_substrates.check_pattern_overlap(
            list(substrate),
            overlap_threshold=overlap_threshold,
        )

    # ── Parse substrate terms ─────────────────────────────────────────────────
    custom_terms = None
    if substrate_terms:
        custom_terms = [t.strip() for t in substrate_terms.split(',')]

    # ── Validate substrates ───────────────────────────────────────────────────
    _section("Validating substrates")
    substrate_family_map = {}
    for sub in substrate:
        try:
            families = parse_substrates.validate_substrate(
                sub, fam_sub_map,
                substrate_terms=custom_terms,
            )
            substrate_family_map[sub] = families

            # Auto-derive patterns if not in patterns file
            existing = pd.read_csv(_PATTERNS_FILE, sep='\t')
            if sub not in existing['substrate'].values:
                click.echo(
                    f"  '{sub}' not in activity_patterns.tsv — "
                    f"deriving patterns automatically..."
                )
                patterns = parse_substrates.auto_derive_patterns(
                    sub, fam_sub_map,
                    substrate_terms=custom_terms,
                )
                parse_substrates.update_patterns_file(
                    sub, patterns, _PATTERNS_FILE,
                    source='auto_derived', reviewed=False,
                )
        except ValueError as e:
            raise click.ClickException(str(e))

    # ── Startup summary ───────────────────────────────────────────────────────
    _section("Pipeline configuration")
    click.echo(f"  Substrates:        {', '.join(substrate)}")
    click.echo(f"  Output:            {output}")
    click.echo(f"  PUL mode:          {pul_mode}")
    click.echo(f"  Min CAZymes/CGC:   {min_substrate_cazymes}")
    click.echo(f"  Threads:           {threads}")
    click.echo(f"  Skip tree:         {skip_tree}")
    click.echo(f"  Skip clinker:      {skip_clinker}")
    click.echo(f"  Reference seqs:    "
               f"{'yes' if ref_seqs else 'no'}")
    click.echo(f"  Reference metadata:"
               f" {'yes' if ref_metadata else 'no'}")

    # ── Per-substrate processing ──────────────────────────────────────────────
    n_total   = len(substrate)
    results   = {}
    log_dir   = os.path.join(output, 'logs')
    os.makedirs(log_dir, exist_ok=True)

    for sub_idx, sub in enumerate(substrate, 1):

        _section(
            f"Substrate: {sub} ({sub_idx}/{n_total})")

        sub_output_dir = os.path.join(output, sub)
        os.makedirs(sub_output_dir, exist_ok=True)

        families      = substrate_family_map[sub]
        activity_file = os.path.join(
            sub_output_dir, f'{sub}_activity_annotated.tsv')
        seq_dir       = os.path.join(sub_output_dir, 'sequences')
        skipped_families = []

        try:

            # ── 1. Classification ─────────────────────────────────────────
            n_steps = 7 if not skip_tree else 5
            if not skip_clinker:
                n_steps += 1
            _step(1, n_steps, "PUL classification")

            substrate_hits, family_hits, overview_df = (
                classify_pul.process_samples(
                    cgc_output_dir=cgc_output_dir,
                    substrate=sub,
                    pul_mode=pul_mode,
                    min_cazymes=min_substrate_cazymes,
                )
            )

            if family_hits.empty:
                _warn(f"No family hits found for {sub} — skipping")
                results[sub] = 'SKIPPED — no hits'
                continue

            n_pul = (family_hits['localisation'] == 'canonical_PUL'
                     ).sum()
            _success(
                f"{len(family_hits)} family hits, "
                f"{n_pul} in canonical PULs"
            )

            # Save classification outputs
            _save_dataframe(
                family_hits,
                os.path.join(sub_output_dir,
                             f'{sub}_family_hits.tsv'),
                'family hits'
            )
            if not substrate_hits.empty:
                _save_dataframe(
                    substrate_hits,
                    os.path.join(sub_output_dir,
                                 f'{sub}_substrate_hits.tsv'),
                    'substrate hits'
                )

            # ── 2. Activity annotation ────────────────────────────────────
            _step(2, n_steps, "Activity annotation")

            family_hits = activity.annotate_hits(
                hits_df=family_hits,
                expasy_file=expasy,
                fam_sub_map=fam_sub_map,
            )

            ref_rows = pd.DataFrame()
            if ref_metadata:
                ref_rows = activity.annotate_references(
                    ref_metadata=ref_metadata,
                    substrate=sub,
                )
                if not ref_rows.empty:
                    _success(
                        f"Appended {len(ref_rows)} reference sequences"
                    )

            all_hits = pd.concat(
                [family_hits, ref_rows], ignore_index=True
            ) if not ref_rows.empty else family_hits

            _save_dataframe(all_hits, activity_file,
                            'activity annotations')

            # Write pattern review report for unreviewed substrates
            existing_patterns = pd.read_csv(_PATTERNS_FILE, sep='\t')
            sub_patterns = existing_patterns[
                existing_patterns['substrate'] == sub
            ]
            if not sub_patterns.empty and not bool(
                    sub_patterns['reviewed'].all()):
                parse_substrates.write_pattern_review_report(
                    substrate=sub,
                    hits_df=family_hits,
                    activity_file=activity_file,
                    patterns=sub_patterns['pattern'].tolist(),
                    output_dir=sub_output_dir,
                )

            # ── 3. Sequence extraction ────────────────────────────────────
            _step(3, n_steps, "Sequence extraction")

            fasta_paths = extract_seqs.extract_sequences(
                cgc_output_dir=cgc_output_dir,
                hits_df=family_hits,
                output_dir=sub_output_dir,
                substrate=sub,
                ref_metadata=ref_metadata,
                ref_seq_dir=ref_seqs,
            )

            if not fasta_paths:
                _warn("No sequences extracted — skipping downstream steps")
                results[sub] = 'SKIPPED — no sequences extracted'
                continue

            # ── 4-6. Placement or de novo tree building ───────────────
            if not skip_tree:
                tree_step = 4
                _step(tree_step, n_steps,
                      "Sequence placement / tree building")

                placement_log_dir = os.path.join(log_dir, sub)
                os.makedirs(placement_log_dir, exist_ok=True)

                # Use de novo for all families if --skip_placement
                # or if epa-ng/hmmer not available
                use_placement = not skip_placement
                if use_placement:
                    try:
                        place_sequences.check_epang()
                        place_sequences.check_hmmer()
                    except place_sequences.ToolNotFoundError as e:
                        _warn(f"EPA placement unavailable ({e}), "
                              f"falling back to de novo for all families")
                        use_placement = False

                if use_placement:
                    tree_results = place_all_families(
                        fasta_paths=fasta_paths,
                        output_dir=sub_output_dir,
                        threads=threads,
                        log_dir=placement_log_dir,
                    )
                    for family, res in tree_results.items():
                        if res['method'] is None:
                            skipped_families.append(
                                (sub, family, 'no sequences'))
                else:
                    # De novo for all families
                    align_dir   = os.path.join(sub_output_dir,
                                               'alignments')
                    trimmed_dir = os.path.join(sub_output_dir,
                                               'trimmed')
                    tree_dir    = os.path.join(sub_output_dir,
                                               'trees')
                    os.makedirs(align_dir,   exist_ok=True)
                    os.makedirs(trimmed_dir, exist_ok=True)
                    os.makedirs(tree_dir,    exist_ok=True)

                    for family, fasta_path in sorted(
                            fasta_paths.items()):
                        try:
                            aligned = align.align(
                                fasta_path=fasta_path,
                                output_path=os.path.join(
                                    align_dir, f'{family}.aln'),
                                threads=threads,
                                log_path=os.path.join(
                                    placement_log_dir,
                                    f'{family}_mafft.log'),
                            )
                            trimmed = trim.trim(
                                alignment_path=aligned,
                                output_path=os.path.join(
                                    trimmed_dir, f'{family}.trim'),
                                log_path=os.path.join(
                                    placement_log_dir,
                                    f'{family}_trimal.log'),
                            )
                            tree.build_tree(
                                trimmed_path=trimmed,
                                output_prefix=os.path.join(
                                    tree_dir, family),
                                threads=threads,
                                log_path=os.path.join(
                                    placement_log_dir,
                                    f'{family}_iqtree.log'),
                            )
                            _success(f"{family} tree built")
                        except TooFewSequencesError as e:
                            _warn(f"Skipping {family}: {e}")
                            skipped_families.append(
                                (sub, family, str(e)))
                        except Exception as e:
                            _warn(f"Failed {family}: {e}")
                            skipped_families.append(
                                (sub, family, str(e)))

            # ── GenBank step number depends on skip_tree ──────────────────
            gbk_step = 4 if skip_tree else 7

            # ── GenBank ───────────────────────────────────────────────────
            _step(gbk_step, n_steps, "GenBank file generation")

            genbank.make_genbank_files(
                cgc_output_dir=cgc_output_dir,
                hits_df=family_hits,
                genomes_dir=genomes or cgc_output_dir,
                output_dir=sub_output_dir,
                substrate=sub,
                tcdb_file=tcdb,
                activity_file=activity_file,
                min_cazymes=min_substrate_cazymes,
                patterns_file=_PATTERNS_FILE,
            )

            # ── iTOL annotations ──────────────────────────────────────────
            itol_step = gbk_step + 1
            _step(itol_step, n_steps,
                  "iTOL annotation files")

            itol.generate_itol_annotations(
                seq_dir=seq_dir,
                output_dir=sub_output_dir,
                substrate=sub,
                colours_file=_COLOURS_FILE,
                activity_file=activity_file,
                ref_metadata=ref_metadata,
                sample_metadata=None,
            )

            # ── Clinker ───────────────────────────────────────────────────
            if not skip_clinker:
                clinker_step = itol_step + 1
                _step(clinker_step, n_steps,
                      "Clinker synteny plot")

                gbk_dir = os.path.join(sub_output_dir, 'genbank')
                gf_file, cm_file = clinker.generate_clinker_inputs(
                    gbk_dir=gbk_dir,
                    output_dir=sub_output_dir,
                    substrate=sub,
                    colours_file=_COLOURS_FILE,
                )

                if gf_file:
                    clinker.run_clinker(
                        gbk_dir=gbk_dir,
                        output_dir=sub_output_dir,
                        substrate=sub,
                        gene_functions_file=gf_file,
                        colour_map_file=cm_file,
                        jobs=threads,
                    )

            results[sub] = 'SUCCESS'

        except click.ClickException:
            raise
        except Exception as e:
            _error(f"Pipeline failed for {sub}: {e}")
            log_path = os.path.join(log_dir, f'{sub}_error.log')
            with open(log_path, 'w') as f:
                traceback.print_exc(file=f)
            _error(f"Full traceback written to {log_path}")
            results[sub] = f'FAILED — see logs/{sub}_error.log'

    # ── Final summary ─────────────────────────────────────────────────────────
    _section("Pipeline complete")

    for sub, status in results.items():
        icon = '✓' if status == 'SUCCESS' else '✗'
        click.echo(f"  {icon} {sub:<20} {status}")

    if skipped_families:
        click.echo(f"\nSkipped families (too few sequences):")
        for sub, family, reason in skipped_families:
            click.echo(f"  {sub}/{family}: {reason}")

    click.echo(f"\nOutputs: {output}")


# ── annotate subcommand ───────────────────────────────────────────────────────

@main.command()
@click.option('--genomes', required=True, type=click.Path(),
              help='Path to directory containing genome FASTA files')
@click.option('--db_dir', required=True, type=click.Path(),
              help='Path to dbCAN database directory')
@click.option('--output', required=True, type=click.Path(),
              help='Output directory')
@click.option('--threads', default=8, show_default=True,
              help='Number of threads')
@click.option('--force', is_flag=True, default=False,
              help='Overwrite existing output files')
def annotate(genomes, db_dir, output, threads, force):
    """Run dbCAN annotation only."""
    _validate_paths(genomes=genomes, db_dir=db_dir)
    cgc_output_dir = run_dbcan.annotate_genomes(
        input_dir=genomes,
        output_dir=output,
        db_dir=db_dir,
        threads=threads,
        force=force,
    )
    _success(f"Annotation complete: {cgc_output_dir}")


# ── classify subcommand ───────────────────────────────────────────────────────

@main.command()
@click.option('--substrate', required=True, multiple=True,
              help='Target substrate(s)')
@click.option('--dbcan_output', required=True, type=click.Path(),
              help='Path to dbCAN cgc_output directory')
@click.option('--db_dir', required=True, type=click.Path(),
              help='Path to dbCAN database directory')
@click.option('--output', required=True, type=click.Path(),
              help='Base output directory')
@click.option('--pul_mode', default='bacteroidetes', show_default=True,
              type=click.Choice(['bacteroidetes', 'generic',
                                 'cazyme_only']),
              help='PUL classification mode')
@click.option('--min_substrate_cazymes', default=2, show_default=True,
              help='Minimum substrate CAZymes required per CGC')
def classify(substrate, dbcan_output, db_dir, output, pul_mode,
             min_substrate_cazymes):
    """Run PUL classification only."""
    _validate_paths(dbcan_output=dbcan_output, db_dir=db_dir)
    fam_sub_map = _get_fam_sub_map(db_dir)

    for sub in substrate:
        sub_output_dir = os.path.join(output, sub)
        os.makedirs(sub_output_dir, exist_ok=True)

        parse_substrates.validate_substrate(sub, fam_sub_map)

        substrate_hits, family_hits, _ = classify_pul.process_samples(
            cgc_output_dir=dbcan_output,
            substrate=sub,
            pul_mode=pul_mode,
            min_cazymes=min_substrate_cazymes,
        )

        if not family_hits.empty:
            _save_dataframe(
                family_hits,
                os.path.join(sub_output_dir,
                             f'{sub}_family_hits.tsv'),
                f'{sub} family hits'
            )
        if not substrate_hits.empty:
            _save_dataframe(
                substrate_hits,
                os.path.join(sub_output_dir,
                             f'{sub}_substrate_hits.tsv'),
                f'{sub} substrate hits'
            )


# ── tree subcommand ───────────────────────────────────────────────────────────

@main.command()
@click.option('--substrate', required=True, multiple=True,
              help='Target substrate(s)')
@click.option('--output', required=True, type=click.Path(),
              help='Base output directory (must contain existing '
                   'sequences/ from a previous run)')
@click.option('--threads', default=8, show_default=True,
              help='Number of threads')
def tree_cmd(substrate, output, threads):
    """Run alignment, trimming and tree building only."""
    for sub in substrate:
        sub_output_dir = os.path.join(output, sub)
        seq_dir        = os.path.join(sub_output_dir, 'sequences')

        if not os.path.exists(seq_dir):
            raise click.ClickException(
                f"No sequences directory found at {seq_dir}. "
                f"Run the full pipeline first."
            )

        tree_dir    = os.path.join(sub_output_dir, 'trees')
        align_dir   = os.path.join(sub_output_dir, 'alignments')
        trimmed_dir = os.path.join(sub_output_dir, 'trimmed')
        log_dir     = os.path.join(output, 'logs')
        os.makedirs(tree_dir,    exist_ok=True)
        os.makedirs(align_dir,   exist_ok=True)
        os.makedirs(trimmed_dir, exist_ok=True)
        os.makedirs(log_dir,     exist_ok=True)

        skipped = []
        for faa_file in sorted(os.listdir(seq_dir)):
            if not faa_file.endswith('.faa'):
                continue

            family      = faa_file.replace('.faa', '')
            fasta_path  = os.path.join(seq_dir, faa_file)
            aligned     = os.path.join(align_dir,
                                       f'{sub}_{family}.aln')
            trimmed     = os.path.join(trimmed_dir,
                                       f'{sub}_{family}.trim')
            tree_prefix = os.path.join(tree_dir,
                                       f'{sub}_{family}')

            try:
                click.echo(f"  {family}...")
                align.align(fasta_path, aligned, threads,
                            os.path.join(log_dir,
                                         f'{sub}_mafft.log'))
                trim.trim(aligned, trimmed,
                          os.path.join(log_dir,
                                       f'{sub}_trimal.log'))
                tree.build_tree(trimmed, tree_prefix, threads,
                                log_path=os.path.join(
                                    log_dir, f'{sub}_iqtree.log'))
                _success(f"{family} done")
            except TooFewSequencesError as e:
                _warn(f"Skipping {family}: {e}")
                skipped.append(family)

        if skipped:
            click.echo(f"\nSkipped: {', '.join(skipped)}")


# ── visualise subcommand ──────────────────────────────────────────────────────

@main.command()
@click.option('--substrate', required=True, multiple=True,
              help='Target substrate(s)')
@click.option('--output', required=True, type=click.Path(),
              help='Base output directory')
@click.option('--ref_metadata', type=click.Path(), default=None,
              help='Path to reference sequence metadata TSV')
@click.option('--colour_config', type=click.Path(), default=None,
              help='Path to existing colour config TSV to use instead '
                   'of regenerating colours')
def visualise(substrate, output, ref_metadata, colour_config):
    """Generate iTOL annotation files only."""
    for sub in substrate:
        sub_output_dir = os.path.join(output, sub)
        seq_dir        = os.path.join(sub_output_dir, 'sequences')
        activity_file  = os.path.join(
            sub_output_dir, f'{sub}_activity_annotated.tsv')

        if not os.path.exists(seq_dir):
            raise click.ClickException(
                f"No sequences directory at {seq_dir}. "
                f"Run the full pipeline first."
            )

        itol.generate_itol_annotations(
            seq_dir=seq_dir,
            output_dir=sub_output_dir,
            substrate=sub,
            colours_file=_COLOURS_FILE,
            activity_file=activity_file,
            ref_metadata=ref_metadata,
            colour_config_path=colour_config,
        )
        _success(f"iTOL annotations written for {sub}")


# ── synteny subcommand ────────────────────────────────────────────────────────

@main.command()
@click.option('--substrate', required=True, multiple=True,
              help='Target substrate(s)')
@click.option('--output', required=True, type=click.Path(),
              help='Base output directory')
@click.option('--threads', default=8, show_default=True,
              help='Number of parallel jobs for clinker')
@click.option('--identity', default=0.3, show_default=True,
              help='Minimum identity threshold for clinker')
def synteny(substrate, output, threads, identity):
    """Run clinker synteny plot only."""
    for sub in substrate:
        sub_output_dir = os.path.join(output, sub)
        gbk_dir        = os.path.join(sub_output_dir, 'genbank')

        if not os.path.exists(gbk_dir):
            raise click.ClickException(
                f"No genbank directory at {gbk_dir}. "
                f"Run the full pipeline first."
            )

        gf_file, cm_file = clinker.generate_clinker_inputs(
            gbk_dir=gbk_dir,
            output_dir=sub_output_dir,
            substrate=sub,
            colours_file=_COLOURS_FILE,
        )

        if gf_file:
            clinker.run_clinker(
                gbk_dir=gbk_dir,
                output_dir=sub_output_dir,
                substrate=sub,
                gene_functions_file=gf_file,
                colour_map_file=cm_file,
                jobs=threads,
                identity=identity,
            )
            _success(f"Clinker plot written for {sub}")


# ── survey subcommand ─────────────────────────────────────────────────────────

@main.command()
@click.option('--dbcan_output', required=True, type=click.Path(),
              help='Path to dbCAN cgc_output directory')
@click.option('--db_dir', required=True, type=click.Path(),
              help='Path to dbCAN database directory')
@click.option('--output', type=click.Path(), default=None,
              help='Optional path to save survey_results.tsv')
@click.option('--pul_mode', default='bacteroidetes', show_default=True,
              type=click.Choice(['bacteroidetes', 'generic',
                                 'cazyme_only']),
              help='PUL classification mode')
@click.option('--min_substrate_cazymes', default=2, show_default=True,
              help='Minimum substrate CAZymes required per CGC')
def survey(dbcan_output, db_dir, output, pul_mode,
           min_substrate_cazymes):
    """
    Survey dbCAN output for all substrates and report hit counts.

    Useful for exploring which substrates are present in your dataset
    before committing to a full pipeline run.
    """
    _validate_paths(dbcan_output=dbcan_output, db_dir=db_dir)
    fam_sub_map = _get_fam_sub_map(db_dir)

    survey_df = parse_substrates.survey_substrates(
        cgc_output_dir=dbcan_output,
        fam_sub_map=fam_sub_map,
        pul_mode=pul_mode,
        min_cazymes=min_substrate_cazymes,
    )

    df_all = pd.read_csv(fam_sub_map, sep='\t', dtype=str)
    df_all.columns = [c.strip() for c in df_all.columns]
    total = df_all['Substrate_curated'].dropna().nunique()

    parse_substrates.print_survey_results(survey_df, total)

    if output and not survey_df.empty:
        os.makedirs(output, exist_ok=True)
        path = os.path.join(output, 'survey_results.tsv')
        survey_df.to_csv(path, sep='\t', index=False)
        _success(f"Survey results saved to {path}")


# ── list-substrates subcommand ────────────────────────────────────────────────

@main.command('list-substrates')
def list_substrates():
    """List all built-in substrates and their CAZyme families."""
    parse_substrates.list_builtin_substrates()




# ── build-reference-db subcommand ─────────────────────────────────────────────

@main.command('build-reference-db')
@click.option('--email', required=True,
              help='Email address for NCBI Entrez API')
@click.option('--api_key', default=None,
              help='NCBI API key (enables 10 requests/second)')
@click.option('--output',
              default='substrate/data/reference_seqs',
              show_default=True,
              help='Output directory for reference sequences')
@click.option('--substrates', multiple=True, default=None,
              help='Substrates to fetch (default: all built-in)')
@click.option('--families', multiple=True, default=None,
              help='Families to fetch (default: all per substrate)')
@click.option('--force', is_flag=True, default=False,
              help='Re-fetch families that already have output files')
def build_reference_db(email, api_key, output, substrates,
                       families, force):
    """
    Build CAZy reference sequence database for EPA placement.

    Downloads characterised enzyme sequences from CAZy and NCBI
    for all built-in substrates. Run this once before building
    reference trees.

    Requires internet access and a valid NCBI email address.
    An API key is recommended for faster downloads (10 req/s vs 3).
    """
    import sys
    sys.path.insert(0, os.path.join(
        os.path.dirname(__file__), '..', 'scripts'))
    from build_reference_db import build_reference_db as _build

    _build(
        email=email,
        api_key=api_key,
        output_dir=output,
        substrates=list(substrates) or None,
        families=list(families) or None,
        force=force,
    )


# ── build-reference-trees subcommand ──────────────────────────────────────────

@main.command('build-reference-trees')
@click.option('--families', multiple=True, default=None,
              help='Families to build trees for '
                   '(default: all in reference_seqs/by_family/)')
@click.option('--max_seqs', type=int, default=None,
              help='Maximum sequences per family tree after subsampling')
@click.option('--threads', default=8, show_default=True,
              help='Number of threads for MAFFT and IQ-TREE2')
@click.option('--min_seqs', type=int, default=4, show_default=True,
              help='Minimum sequences required to build a reference tree')
@click.option('--force', is_flag=True, default=False,
              help='Rebuild trees that already exist')
@click.option('--ref_seqs_dir',
              default='substrate/data/reference_seqs/by_family',
              show_default=True,
              help='Path to by_family/ reference sequences directory')
@click.option('--output',
              default='substrate/data/reference_trees',
              show_default=True,
              help='Output directory for reference trees')
def build_reference_trees(families, max_seqs, threads, min_seqs,
                          force, ref_seqs_dir, output):
    """
    Build IQ-TREE2 reference trees from CAZy characterised sequences.

    Run this after build-reference-db. Trees are stored in
    substrate/data/reference_trees/ and used automatically
    during pipeline runs for EPA sequence placement.

    Families with fewer than --min_seqs sequences are skipped
    and will use de novo tree building during pipeline runs.
    """
    import sys
    sys.path.insert(0, os.path.join(
        os.path.dirname(__file__), '..', 'scripts'))
    from build_reference_trees import main as _main
    import sys as _sys

    # Build argv for the script's argparse
    argv = [
        '--threads',  str(threads),
        '--min_seqs', str(min_seqs),
        '--ref_seqs_dir', ref_seqs_dir,
        '--output',   output,
    ]
    if families:
        argv += ['--families'] + list(families)
    if max_seqs:
        argv += ['--max_seqs_per_family', str(max_seqs)]
    if force:
        argv.append('--force')

    _sys.argv = ['build_reference_trees'] + argv
    _main()
