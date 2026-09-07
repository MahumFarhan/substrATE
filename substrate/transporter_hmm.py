"""
HMM-based SusC/SusD transporter detection, supplementing the existing
TCDB-based check in classify_pul.py.

Background: TCDB's SusD coverage is sparse (13 reference sequences vs 118
for SusC), so TCDB-only classification under-calls canonical_PUL for CGCs
whose transporter gene is a real SusD but wasn't caught by TCDB's diamond
search. This module runs two additional, well-established HMM profiles
directly against each genome's already-predicted protein set
(uniInput.faa, produced by dbCAN's own CGC-Finder step — no
re-annotation or genome re-processing needed):

  - TIGR04056 (TIGRFAM) — SusC-like TonB-dependent transporter / porin
  - PF07980   (Pfam)    — SusD/RagB family

Note: TIGR04056 is the SusC-side model, not SusD — there is no
established TIGRFAM model specifically for SusD; PF07980 is the
standard model for that side. See classify_pul.py's docstring for how
these combine with the existing TCDB check.
"""

import os
import subprocess
import pandas as pd


SUSC_HMM_ACCESSION = 'TIGR04056'
SUSD_HMM_ACCESSION = 'PF07980'


def run_hmmsearch(hmm_path, faa_path, tblout_path, cut_tc=True):
    """
    Run hmmsearch of a single HMM profile against a protein FASTA,
    using the model's built-in trusted cutoff (gathering threshold) by
    default rather than a generic e-value — the correct convention for
    both TIGRFAM and Pfam models.

    Args:
        hmm_path:    path to a single-model .hmm file (e.g. TIGR04056.hmm)
        faa_path:    path to the protein FASTA to search (uniInput.faa)
        tblout_path: path to write hmmsearch's --tblout table to
        cut_tc:      use --cut_tc (trusted cutoff) if True; if the HMM
                     file lacks TC lines this will cause hmmsearch to
                     error, so callers should verify the profile source
                     includes cutoffs (both TIGRFAM and Pfam-A do).

    Raises:
        subprocess.CalledProcessError if hmmsearch exits non-zero.
    """
    cmd = ['hmmsearch']
    if cut_tc:
        cmd.append('--cut_tc')
    cmd += ['--tblout', tblout_path, '--noali', hmm_path, faa_path]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def parse_tblout_gene_ids(tblout_path):
    """
    Parse an hmmsearch --tblout file and return the set of target
    (query protein) names that produced a hit — i.e. every gene ID that
    passed the model's cutoff.

    --tblout is whitespace-delimited with '#'-prefixed comment/header
    lines; the target name is the first field of each data line.
    """
    gene_ids = set()
    if not os.path.exists(tblout_path):
        return gene_ids
    with open(tblout_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            gene_ids.add(line.split()[0])
    return gene_ids


def annotate_sample(sample, cgc_output_dir, susc_hmm_path, susd_hmm_path,
                    log_dir=None):
    """
    Run both SusC and SusD HMM searches for one sample and return a
    small DataFrame of hits with provenance.

    Args:
        sample:          sample name (matches the 'output_{sample}'
                         directory name under cgc_output_dir)
        cgc_output_dir:  path to the dbCAN cgc_output directory
        susc_hmm_path:   path to TIGR04056.hmm
        susd_hmm_path:   path to PF07980.hmm
        log_dir:         optional directory to write hmmsearch tblout
                         files to for debugging; defaults to a temp
                         location inside the sample's own output dir

    Returns:
        DataFrame with columns: sample, Protein ID, transporter_role,
        hmm_accession, hmm_source. Empty DataFrame (correct columns,
        zero rows) if uniInput.faa is missing for this sample.
    """
    sample_dir = os.path.join(cgc_output_dir, f'output_{sample}')
    faa_path = os.path.join(sample_dir, 'uniInput.faa')

    columns = ['sample', 'Protein ID', 'transporter_role',
              'hmm_accession', 'hmm_source']

    if not os.path.exists(faa_path):
        return pd.DataFrame(columns=columns)

    out_dir = log_dir or sample_dir
    os.makedirs(out_dir, exist_ok=True)

    rows = []

    susc_tblout = os.path.join(out_dir, f'{sample}_TIGR04056.tblout')
    run_hmmsearch(susc_hmm_path, faa_path, susc_tblout, cut_tc=True)
    for gene_id in parse_tblout_gene_ids(susc_tblout):
        rows.append({
            'sample': sample, 'Protein ID': gene_id,
            'transporter_role': 'SusC',
            'hmm_accession': SUSC_HMM_ACCESSION,
            'hmm_source': 'TIGRFAM',
        })

    susd_tblout = os.path.join(out_dir, f'{sample}_PF07980.tblout')
    run_hmmsearch(susd_hmm_path, faa_path, susd_tblout, cut_tc=True)
    for gene_id in parse_tblout_gene_ids(susd_tblout):
        rows.append({
            'sample': sample, 'Protein ID': gene_id,
            'transporter_role': 'SusD',
            'hmm_accession': SUSD_HMM_ACCESSION,
            'hmm_source': 'Pfam',
        })

    return pd.DataFrame(rows, columns=columns)


def annotate_all_samples(cgc_output_dir, susc_hmm_path, susd_hmm_path):
    """
    Run HMM-based transporter detection for every sample directory found
    in cgc_output_dir.

    Returns:
        Combined DataFrame across all samples (same columns as
        annotate_sample). Write this to a single
        'transporter_hmm_hits.tsv' at the cgc_output_dir level, or
        wherever classify_pul.py's caller expects to load it from —
        this is the one piece of wiring I haven't finalized yet (see
        chat) since I don't yet know the exact path convention
        process_samples()/cli.py should read it from.
    """
    all_hits = []
    for entry in sorted(os.listdir(cgc_output_dir)):
        if not entry.startswith('output_'):
            continue
        sample = entry[len('output_'):]
        hits = annotate_sample(sample, cgc_output_dir,
                               susc_hmm_path, susd_hmm_path)
        if not hits.empty:
            all_hits.append(hits)

    columns = ['sample', 'Protein ID', 'transporter_role',
              'hmm_accession', 'hmm_source']
    if not all_hits:
        return pd.DataFrame(columns=columns)
    return pd.concat(all_hits, ignore_index=True)
