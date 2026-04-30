"""
Build reference sequence database from CAZy characterised enzymes.

For each substrate in FAMILY_MAP, fetches characterised enzyme
sequences from CAZy and NCBI, organised into per-substrate,
per-family FASTA files for use in EPA-RAxML reference trees.

Output structure:
    substrate/data/reference_seqs/
    ├── laminarin/
    │   ├── GH16.faa
    │   ├── GH17.faa
    │   └── GH55.faa
    ├── alginate/
    │   └── PL7.faa
    └── reference_metadata.tsv

Usage:
    python scripts/build_reference_db.py \\
        --email mahum.farhan@gmail.com \\
        --api_key 0ea2c3495834b7a2844a6faa897166179308 \\
        --output substrate/data/reference_seqs \\
        [--substrates laminarin alginate] \\
        [--families GH16 GH17] \\
        [--force]

If --substrates is omitted, all built-in substrates are processed.
If --families is omitted, all families for each substrate are fetched.
Use --force to re-fetch families that already have output files.
"""
import os
import re
import time
import json
import argparse
import requests
import pandas as pd
from io import StringIO
from bs4 import BeautifulSoup
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# Add parent directory to path so we can import substrate modules
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))

from substrate.classify_pul import FAMILY_MAP


# ── Constants ─────────────────────────────────────────────────────────────────

CAZY_BASE_URL   = 'https://www.cazy.org'
NCBI_BATCH_SIZE = 200    # Max accessions per Entrez fetch
CAZY_SLEEP      = 2.0    # Seconds between CAZy requests (be polite)
NCBI_SLEEP_AUTH = 0.15   # Seconds between NCBI requests with API key
NCBI_SLEEP_ANON = 0.4    # Seconds between NCBI requests without key

# Cache directory for raw CAZy HTML and NCBI responses
# Allows interrupted runs to resume without re-fetching
CACHE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '.cazy_cache')


# ── CAZy HTML parsing ─────────────────────────────────────────────────────────

def fetch_cazy_characterised(family, session, cache=True):
    """
    Fetch and parse the characterised enzyme table for a CAZy family.

    Args:
        family:  CAZy family string (e.g. 'GH16')
        session: requests.Session object
        cache:   if True, cache raw HTML to avoid re-fetching

    Returns:
        list of dicts with keys: accessions, organism, ec_numbers,
        subfamily, protein_name
        Returns empty list if page not found or no table present.
    """
    url        = f"{CAZY_BASE_URL}/{family}_characterized.html"
    cache_path = os.path.join(CACHE_DIR, f"{family}_char.html")

    # Try cache first
    if cache and os.path.exists(cache_path):
        print(f"    Using cached page for {family}")
        with open(cache_path, 'r', encoding='utf-8') as f:
            html = f.read()
    else:
        print(f"    Fetching {url}")
        try:
            response = session.get(url, timeout=30)
            if response.status_code == 404:
                print(f"    WARNING: {family} not found on CAZy")
                return []
            response.raise_for_status()
            html = response.text
            time.sleep(CAZY_SLEEP)

            if cache:
                os.makedirs(CACHE_DIR, exist_ok=True)
                with open(cache_path, 'w', encoding='utf-8') as f:
                    f.write(html)
        except requests.RequestException as e:
            print(f"    WARNING: Could not fetch {family}: {e}")
            return []

    return _parse_cazy_table(html, family)


def _parse_cazy_table(html, family):
    """
    Parse a CAZy characterised enzyme HTML table.

    The table structure on CAZy pages:
      Protein Name | EC# | Reference | Organism | GenBank |
      UniProt | PDB/3D | Subf

    Accessions are extracted from anchor tags in the GenBank cell
    rather than regex on text, to correctly capture 9-digit WP_
    accessions and other edge cases.

    Args:
        html:   raw HTML string
        family: family name for logging

    Returns:
        list of dicts with keys: accessions, organism, ec_numbers,
        subfamily, protein_name
    """
    soup  = BeautifulSoup(html, 'html.parser')
    rows  = []

    # CAZy tables have class 'listing'
    tables = soup.find_all('table', {'class': 'listing'})
    if not tables:
        print(f"    WARNING: No characterised table found for {family}")
        return []

    # Valid NCBI protein accession patterns:
    # WP_123456789.1  (9 digits, RefSeq non-redundant)
    # NP_123456.1     (6 digits, RefSeq)
    # XP_001234567.1  (9 digits, predicted)
    # AAC25554.2      (2 letters + 5-8 digits)
    # ABN51452.1      (3 letters + 5-8 digits)
    ACC_PATTERN = re.compile(
        r'(?:WP|NP|XP|YP|AP|ZP)_[0-9]{6,9}(?:\.[0-9]+)?'
        r'|[A-Z]{2,3}[0-9]{5,8}(?:\.[0-9]+)?'
    )

    for table in tables:
        for tr in table.find_all('tr'):
            cells = tr.find_all('td')
            if len(cells) < 5:
                continue

            protein_name = cells[0].get_text(separator=' ',
                                              strip=True)
            ec_text      = cells[1].get_text(separator=' ',
                                             strip=True)
            organism     = cells[3].get_text(separator=' ',
                                             strip=True)
            genbank_cell = cells[4]
            subfamily    = cells[-1].get_text(strip=True) \
                if len(cells) >= 8 else ''

            # Extract accessions from anchor tags first (most reliable)
            accessions = []
            for a in genbank_cell.find_all('a'):
                href = a.get('href', '')
                text = a.get_text(strip=True)
                # Try href first (contains full accession)
                for match in ACC_PATTERN.findall(href):
                    if match not in accessions:
                        accessions.append(match)
                # Also try link text
                for match in ACC_PATTERN.findall(text):
                    if match not in accessions:
                        accessions.append(match)

            # Fallback: regex on cell text if no anchors found
            if not accessions:
                cell_text = genbank_cell.get_text(
                    separator=' ', strip=True)
                accessions = ACC_PATTERN.findall(cell_text)

            if not accessions:
                continue

            # Parse EC numbers
            ec_numbers = re.findall(
                r'\d+\.\d+\.\d+\.[\d\-]+', ec_text)

            rows.append({
                'accessions':   accessions,
                'organism':     organism,
                'ec_numbers':   ec_numbers,
                'subfamily':    subfamily,
                'protein_name': protein_name,
            })

    print(f"    Found {len(rows)} characterised entries for {family}")
    return rows


# ── Sequence subsampling ─────────────────────────────────────────────────────

def subsample_by_subfamily_diversity(records_dict, metadata_rows,
                                     max_seqs):
    """
    Subsample sequences to max_seqs using Option C strategy:
    proportional representation across subfamilies, then by
    organism diversity within each subfamily.

    This ensures all subfamilies are represented even when capping,
    which is important since subfamilies often correspond to distinct
    substrate specificities.

    Args:
        records_dict:  dict mapping accession to SeqRecord
        metadata_rows: list of metadata dicts for these sequences
        max_seqs:      maximum number of sequences to keep

    Returns:
        tuple of (subsampled_records_dict, subsampled_metadata_rows)
    """
    if len(records_dict) <= max_seqs:
        return records_dict, metadata_rows

    # Group metadata rows by subfamily
    subfamily_groups = {}
    for row in metadata_rows:
        subfam = row['subfamily'] or 'unknown'
        if subfam not in subfamily_groups:
            subfamily_groups[subfam] = []
        for acc in row['accessions']:
            if acc in records_dict or                     acc.split('.')[0] in records_dict:
                subfamily_groups[subfam].append(row)
                break

    n_subfamilies = len(subfamily_groups)
    if n_subfamilies == 0:
        return records_dict, metadata_rows

    # Allocate slots proportionally across subfamilies
    # Each subfamily gets at least 1 slot, remainder distributed
    # proportionally by subfamily size
    base_slots  = max(1, max_seqs // n_subfamilies)
    total_sizes = {sf: len(rows)
                   for sf, rows in subfamily_groups.items()}
    total_seqs  = sum(total_sizes.values())

    slots = {}
    for sf, size in total_sizes.items():
        proportion  = size / total_seqs
        slots[sf]   = max(1, int(max_seqs * proportion))

    # Adjust to hit exact max_seqs target
    while sum(slots.values()) > max_seqs:
        largest = max(slots, key=lambda s: slots[s])
        slots[largest] -= 1
    while sum(slots.values()) < max_seqs:
        # Give extra slots to largest subfamilies
        for sf in sorted(slots, key=lambda s: -total_sizes[s]):
            if slots[sf] < total_sizes[sf]:
                slots[sf] += 1
                if sum(slots.values()) >= max_seqs:
                    break

    # Within each subfamily, select by organism diversity
    # (one sequence per organism first, then fill remaining slots)
    selected_accs     = set()
    selected_metadata = []

    for subfam, n_slots in slots.items():
        rows       = subfamily_groups[subfam]
        # Group by organism
        by_organism = {}
        for row in rows:
            org = row['organism']
            if org not in by_organism:
                by_organism[org] = []
            by_organism[org].append(row)

        # Round-robin across organisms for diversity
        selected_rows = []
        org_lists     = list(by_organism.values())
        i = 0
        while len(selected_rows) < n_slots and any(org_lists):
            org_idx = i % len(org_lists)
            if org_lists[org_idx]:
                row = org_lists[org_idx].pop(0)
                # Check we have a sequence for this row
                for acc in row['accessions']:
                    base = acc.split('.')[0]
                    if acc in records_dict or base in records_dict:
                        if acc not in selected_accs:
                            selected_accs.add(acc)
                            selected_rows.append(row)
                            break
            i += 1
            # Avoid infinite loop if all lists exhausted
            if all(not ol for ol in org_lists):
                break

        selected_metadata.extend(selected_rows)

    # Build subsampled records dict
    subsampled = {}
    for row in selected_metadata:
        for acc in row['accessions']:
            base = acc.split('.')[0]
            if acc in records_dict:
                subsampled[acc] = records_dict[acc]
                break
            elif base in records_dict:
                subsampled[base] = records_dict[base]
                break

    print(f"    Subsampled {len(records_dict)} -> "
          f"{len(subsampled)} sequences "
          f"({n_subfamilies} subfamilies, "
          f"max {max_seqs} per family)")

    return subsampled, selected_metadata


# ── NCBI sequence fetching ────────────────────────────────────────────────────

def fetch_sequences_ncbi(accessions, sleep_time, cache=True):
    """
    Fetch protein sequences from NCBI for a list of accessions.

    Uses batch fetching (up to NCBI_BATCH_SIZE per request) to
    minimise API calls. Caches results to allow resume on failure.

    Args:
        accessions: list of GenBank accession strings
        sleep_time: seconds to sleep between requests
        cache:      if True, cache fetched sequences

    Returns:
        dict mapping accession to SeqRecord
    """
    results    = {}
    to_fetch   = []
    cache_file = os.path.join(CACHE_DIR, 'ncbi_sequences.faa')

    # Load cached sequences
    cached = {}
    if cache and os.path.exists(cache_file):
        for record in SeqIO.parse(cache_file, 'fasta'):
            cached[record.id.split('.')[0]] = record
            cached[record.id]               = record

    for acc in accessions:
        acc_base = acc.split('.')[0]
        if acc in cached:
            results[acc] = cached[acc]
        elif acc_base in cached:
            results[acc] = cached[acc_base]
        else:
            to_fetch.append(acc)

    if not to_fetch:
        return results

    print(f"    Fetching {len(to_fetch)} sequences from NCBI "
          f"({len(accessions) - len(to_fetch)} cached)...")

    new_records = {}
    for i in range(0, len(to_fetch), NCBI_BATCH_SIZE):
        batch = to_fetch[i:i + NCBI_BATCH_SIZE]
        try:
            handle = Entrez.efetch(
                db='protein',
                id=','.join(batch),
                rettype='fasta',
                retmode='text'
            )
            for record in SeqIO.parse(handle, 'fasta'):
                new_records[record.id] = record
                # Also index by base accession
                base = record.id.split('.')[0]
                new_records[base] = record
            handle.close()
            time.sleep(sleep_time)

            print(f"    Fetched batch {i//NCBI_BATCH_SIZE + 1}/"
                  f"{(len(to_fetch)-1)//NCBI_BATCH_SIZE + 1}")

        except Exception as e:
            print(f"    WARNING: NCBI fetch failed for batch "
                  f"starting at {i}: {e}")
            time.sleep(sleep_time * 3)

    # Append new records to cache
    if cache and new_records:
        os.makedirs(CACHE_DIR, exist_ok=True)
        with open(cache_file, 'a') as f:
            SeqIO.write(list(new_records.values()), f, 'fasta')

    # Map original accessions to records
    for acc in to_fetch:
        acc_base = acc.split('.')[0]
        if acc in new_records:
            results[acc] = new_records[acc]
        elif acc_base in new_records:
            results[acc] = new_records[acc_base]
        else:
            print(f"    WARNING: Could not fetch {acc}")

    return results


# ── Output writing ────────────────────────────────────────────────────────────

def write_family_fasta(records, metadata_rows, output_dir,
                       substrate, family):
    """
    Write sequences for one family to a FASTA file.

    Sequence headers use the format:
        >accession|family|subfamily|organism

    Args:
        records:      dict mapping accession to SeqRecord
        metadata_rows: list of metadata dicts for these sequences
        output_dir:   base reference_seqs directory
        substrate:    substrate name
        family:       family name

    Returns:
        number of sequences written
    """
    sub_dir = os.path.join(output_dir, substrate)
    os.makedirs(sub_dir, exist_ok=True)

    out_path   = os.path.join(sub_dir, f"{family}.faa")
    written    = []
    seen_accs  = set()

    for row in metadata_rows:
        for acc in row['accessions']:
            if acc in seen_accs:
                continue
            acc_base = acc.split('.')[0]

            record = records.get(acc) or records.get(acc_base)
            if record is None:
                continue

            seen_accs.add(acc)
            seen_accs.add(acc_base)

            subfam  = row['subfamily'] or family
            org     = row['organism'].replace(' ', '_')[:40]
            new_rec = SeqRecord(
                record.seq,
                id=f"{acc}|{family}|{subfam}|{org}",
                description=''
            )
            written.append(new_rec)

    if written:
        SeqIO.write(written, out_path, 'fasta')
        print(f"    Wrote {len(written)} sequences -> {out_path}")
    else:
        print(f"    WARNING: No sequences written for {family}")

    return len(written)


# ── Main ──────────────────────────────────────────────────────────────────────

def build_reference_db(email, api_key, output_dir,
                       substrates=None, families=None, force=False):
    """
    Build the reference sequence database for all specified substrates.

    Args:
        email:      NCBI email (required)
        api_key:    NCBI API key (optional but recommended)
        output_dir: path to write reference_seqs/
        substrates: list of substrate names to process (None = all)
        families:   list of family names to process (None = all)
        force:      if True, re-fetch existing files
    """
    Entrez.email   = email
    Entrez.api_key = api_key if api_key else None
    sleep_time     = NCBI_SLEEP_AUTH if api_key else NCBI_SLEEP_ANON

    if api_key:
        print(f"Using NCBI API key — {int(1/sleep_time)} "
              f"requests/second")
    else:
        print(f"No API key — {int(1/sleep_time)} requests/second. "
              f"Consider providing --api_key for faster downloads.")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(CACHE_DIR,  exist_ok=True)

    # Determine which substrates to process
    target_substrates = substrates or list(FAMILY_MAP.keys())
    print(f"\nProcessing {len(target_substrates)} substrates: "
          f"{', '.join(sorted(target_substrates))}\n")

    session       = requests.Session()
    session.headers.update({'User-Agent': 'substrATE/0.1.0'})
    all_metadata  = []
    total_written = 0

    for substrate in sorted(target_substrates):
        if substrate not in FAMILY_MAP:
            print(f"WARNING: '{substrate}' not in FAMILY_MAP, "
                  f"skipping")
            continue

        target_families = families or FAMILY_MAP[substrate]
        print(f"\n{'─'*50}")
        print(f"Substrate: {substrate} "
              f"({len(target_families)} families)")
        print(f"{'─'*50}")

        for family in sorted(target_families):
            out_path = os.path.join(output_dir, substrate,
                                    f"{family}.faa")

            if os.path.exists(out_path) and not force:
                existing = sum(
                    1 for _ in SeqIO.parse(out_path, 'fasta'))
                print(f"  {family}: skipping "
                      f"({existing} seqs exist, use --force "
                      f"to re-fetch)")
                continue

            print(f"\n  {family}:")

            # Fetch CAZy characterised table
            char_rows = fetch_cazy_characterised(
                family, session, cache=True)

            if not char_rows:
                print(f"    No characterised entries found")
                continue

            # Collect all accessions for batch NCBI fetch
            all_accessions = []
            for row in char_rows:
                all_accessions.extend(row['accessions'])
            all_accessions = list(dict.fromkeys(all_accessions))

            # Fetch sequences from NCBI
            records = fetch_sequences_ncbi(
                all_accessions, sleep_time, cache=True)

            # Write FASTA
            n_written = write_family_fasta(
                records, char_rows, output_dir, substrate, family)
            total_written += n_written

            # Collect metadata
            for row in char_rows:
                for acc in row['accessions']:
                    if records.get(acc) or \
                            records.get(acc.split('.')[0]):
                        all_metadata.append({
                            'accession':    acc,
                            'family':       family,
                            'subfamily':    row['subfamily']
                                            or family,
                            'substrate':    substrate,
                            'organism':     row['organism'],
                            'ec_numbers':   ','.join(
                                row['ec_numbers']),
                            'protein_name': row['protein_name'],
                            'label':        (
                                f"{row['organism'][:30]} "
                                f"{row['protein_name'][:30]}"
                            ),
                            'source':       'CAZy_characterised',
                        })

    # Write master metadata TSV
    if all_metadata:
        meta_df   = pd.DataFrame(all_metadata).drop_duplicates(
            subset=['accession'])
        meta_path = os.path.join(output_dir,
                                 'reference_metadata.tsv')
        meta_df.to_csv(meta_path, sep='\t', index=False)
        print(f"\nMetadata written to {meta_path} "
              f"({len(meta_df)} sequences)")

    # ── Merge into global per-family FASTAs ──────────────────────
    print(f"\nMerging into global per-family FASTAs...")
    by_family_dir = os.path.join(output_dir, 'by_family')
    os.makedirs(by_family_dir, exist_ok=True)

    # Collect all sequences grouped by family
    family_records = {}
    for substrate_dir in sorted(os.listdir(output_dir)):
        full_path = os.path.join(output_dir, substrate_dir)
        if not os.path.isdir(full_path) or                 substrate_dir in ('by_family',):
            continue
        for faa_file in sorted(os.listdir(full_path)):
            if not faa_file.endswith('.faa'):
                continue
            family   = faa_file[:-4]
            faa_path = os.path.join(full_path, faa_file)
            if family not in family_records:
                family_records[family] = {}
            for record in SeqIO.parse(faa_path, 'fasta'):
                # Deduplicate by accession (first part of ID)
                acc = record.id.split('|')[0]
                if acc not in family_records[family]:
                    family_records[family][acc] = record

    # Write merged FASTAs (with optional subsampling)
    print(f"\nGlobal family FASTAs:")
    for family, records_dict in sorted(family_records.items()):
        out_path = os.path.join(by_family_dir, f'{family}.faa')

        SeqIO.write(list(records_dict.values()), out_path, 'fasta')
        n = len(records_dict)
        print(f"  {family:<10} {n:>5} sequences -> {out_path}")

    print(f"\n{'='*50}")
    print(f"Done. Total sequences written: {total_written}")
    print(f"Reference sequences: {output_dir}")
    print(f"Global family FASTAs: {by_family_dir}")
    print(f"Cache: {CACHE_DIR}")
    print(f"\nNext step: run scripts/build_reference_trees.py")
    print(f"{'='*50}")


def main():
    parser = argparse.ArgumentParser(
        description='Build CAZy reference sequence database'
    )
    parser.add_argument(
        '--email', required=True,
        help='Email address for NCBI Entrez API (required by NCBI)'
    )
    parser.add_argument(
        '--api_key', default=None,
        help='NCBI API key (optional, enables 10 requests/second)'
    )
    parser.add_argument(
        '--output',
        default='substrate/data/reference_seqs',
        help='Output directory for reference sequences'
    )
    parser.add_argument(
        '--substrates', nargs='+', default=None,
        help='Substrates to process (default: all built-in)'
    )
    parser.add_argument(
        '--families', nargs='+', default=None,
        help='Families to process (default: all for each substrate)'
    )
    parser.add_argument(
        '--force', action='store_true',
        help='Re-fetch families that already have output files'
    )
    args = parser.parse_args()

    build_reference_db(
        email=args.email,
        api_key=args.api_key,
        output_dir=args.output,
        substrates=args.substrates,
        families=args.families,
        force=args.force,
    )


if __name__ == '__main__':
    main()
