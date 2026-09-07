#!/usr/bin/env bash
# Downloads the full TIGRFAM and Pfam-A HMM libraries, extracts just the
# SusC (TIGR04056) and SusD (PF07980) profiles via hmmfetch, and cleans up
# the large source libraries afterward — leaving only the two small HMM
# files substrATE actually needs (--susc_hmm / --susd_hmm).
#
# Requires: HMMER (hmmfetch, hmmpress) already installed and on PATH.
# Run from wherever you want the reference databases to live, e.g.
# ~/substrATE/databases/transporters/

set -euo pipefail

mkdir -p transporters
cd transporters

echo "=== Downloading TIGRFAM release 15.0 (full library) ==="
wget -N https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/TIGRFAMs_15.0_HMM.LIB.gz
gunzip -kf TIGRFAMs_15.0_HMM.LIB.gz

echo ""
echo "=== Downloading Pfam-A (current release, full library) ==="
echo "    This is a large file (~250MB compressed) — one-time download."
wget -N https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip -kf Pfam-A.hmm.gz

echo ""
echo "=== Indexing libraries for hmmfetch ==="
hmmfetch --index TIGRFAMs_15.0_HMM.LIB
hmmfetch --index Pfam-A.hmm

echo ""
echo "=== Extracting SusC (TIGR04056) and SusD (PF07980) ==="
hmmfetch TIGRFAMs_15.0_HMM.LIB TIGR04056 > TIGR04056.hmm
hmmfetch Pfam-A.hmm PF07980 > PF07980.hmm

# Verify both extractions actually got a real model, not an empty file
# (hmmfetch silently writes an empty file if the accession isn't found —
# it does not error out, so this check matters).
for f in TIGR04056.hmm PF07980.hmm; do
    if ! grep -q "^HMM" "$f"; then
        echo "ERROR: $f does not look like a valid HMM — extraction may have failed."
        echo "       Check that the accession still exists in this release."
        exit 1
    fi
done

echo ""
echo "=== Cleaning up large source libraries ==="
echo "    Keeping only TIGR04056.hmm and PF07980.hmm."
rm -f TIGRFAMs_15.0_HMM.LIB TIGRFAMs_15.0_HMM.LIB.gz \
      TIGRFAMs_15.0_HMM.LIB.ssi \
      Pfam-A.hmm Pfam-A.hmm.gz Pfam-A.hmm.ssi

echo ""
echo "Done. Use these with substrate run / substrate classify:"
echo "  --susc_hmm $(pwd)/TIGR04056.hmm"
echo "  --susd_hmm $(pwd)/PF07980.hmm"
echo ""
echo "TIGRFAM (JCVI/NCBI) and Pfam-A (EBI) are both freely redistributable"
echo "under their respective open licenses (CC-BY / CC0) — cite both"
echo "databases in any resulting publication per their usage terms."