# API Integration Specs: UniProt REST + IEDB MHC-II

Researched: 2026-04-02
Tested against live APIs with Python requests from venv.

---

## 1. UniProt REST API -- Protein Sequence Retrieval

### 1.1 Direct Accession Lookup (FASTA)

**Endpoint:** `GET https://rest.uniprot.org/uniprotkb/{accession}.fasta`

**Response:**
- Status 200, Content-Type: `text/plain;format=fasta`
- Standard FASTA: header line starting with `>`, then sequence lines

**Tested accessions:**

| Accession | Actual Protein (check!) | Organism | Seq Length |
|-----------|------------------------|----------|------------|
| Q25327 | Membrane-bound acid phosphatase | L. mexicana (5665) | 516 aa |
| P15711 | 104 kDa microneme/rhoptry antigen | T. parva (5875) | 924 aa |
| Q25223 | DNA-directed RNA polymerase (Fragment) | L. amazonensis (5659) | 422 aa |
| Q25292 | Heat shock 70 kDa protein (Fragment) | L. infantum (5671) | 653 aa |

**IMPORTANT:** The accessions provided in the task description do NOT map to the
expected proteins. P15711 is NOT KMP-11 (it is a Theileria parva protein), Q25327
is NOT an A2 antigen, etc. The pipeline should use search-based lookup or a curated
accession mapping. See section 1.2.

**Error handling:**
- Invalid accession format -> Status 400, body: `"The 'accession' value has invalid format..."`
- Non-existent valid format -> Status 400 (same message)

**Rate limits:**
- No explicit rate-limit headers returned (no X-RateLimit-*, Retry-After, etc.)
- UniProt docs recommend max 1 request per 100ms for automated access
- 5 rapid sequential requests completed in ~4.2s with no throttling observed
- Recommendation: use 0.5s delay between requests as courtesy

### 1.2 Search by Protein Name (FASTA)

**Endpoint:** `GET https://rest.uniprot.org/uniprotkb/search?query={query}&format=fasta&size={n}`

**Query syntax:** Lucene-style field queries. Key fields:
- `protein_name:X` -- search protein description
- `gene:X` -- search gene name
- `organism_id:5671` -- L. infantum taxonomy ID

**Tested searches (organism_id:5671 = L. infantum):**

| Search | Results | Top Hit |
|--------|---------|---------|
| `protein_name:LiHyp1` | 0 | (use `gene:` instead, or broader search) |
| `gene:LiHyp1` | 0 | Not in UniProt by this name |
| `protein_name:HSP70` | 9 | A0A6L0XI79 (Hsp70, 654 aa) |
| `protein_name:A2` | 3+ | A4HZU7 (A2 protein, L. infantum) |
| `gene:LACK` | 3+ | J9XRQ3 (Antigen LACK, L. infantum) |
| `protein_name:GP63` | 0 | (try `gene:GP63` or `protein_name:leishmanolysin`) |
| `protein_name:KMP-11` | 0 | (try `protein_name:kinetoplastid membrane protein`) |

**Pagination:**
- `X-Total-Results` header gives total count
- `Link` header with `rel="next"` provides cursor-based pagination URL
- Default page size: 25; set with `&size=N` (max 500)

### 1.3 Search with JSON format

**Endpoint:** `GET https://rest.uniprot.org/uniprotkb/search?query={query}&format=json&size={n}`

Returns structured JSON with `results` array. Each entry includes:
- `primaryAccession` -- e.g., "A0A6L0XI79"
- `proteinDescription.submissionNames[].fullName.value` -- protein name
- `sequence.value` -- full amino acid sequence
- `sequence.length` -- sequence length
- `organism.taxonId` -- taxonomy ID

### 1.4 Recommended Implementation

```python
import requests
from typing import Optional

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
REQUEST_DELAY = 0.5  # seconds between requests


def fetch_sequence_by_accession(accession: str) -> Optional[str]:
    """Fetch protein sequence from UniProt by accession number.

    Returns the amino acid sequence string, or None if not found.
    """
    url = f"{UNIPROT_BASE}/{accession}.fasta"
    resp = requests.get(url, timeout=30)

    if resp.status_code != 200:
        return None

    lines = resp.text.strip().split("\n")
    # Skip header line(s) starting with '>'
    sequence = "".join(line for line in lines if not line.startswith(">"))
    return sequence


def search_sequences(
    protein_name: str,
    organism_id: int = 5671,
    max_results: int = 5,
) -> list[dict]:
    """Search UniProt for protein sequences by name and organism.

    Args:
        protein_name: Protein name to search (e.g., "HSP70", "A2").
        organism_id: NCBI taxonomy ID (5671 = L. infantum).
        max_results: Maximum entries to return.

    Returns:
        List of dicts with keys: accession, name, sequence, length.
    """
    query = f"protein_name:{protein_name} AND organism_id:{organism_id}"
    params = {"query": query, "format": "json", "size": max_results}
    resp = requests.get(f"{UNIPROT_BASE}/search", params=params, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    results = []
    for entry in data.get("results", []):
        # Extract protein name from various possible locations
        desc = entry.get("proteinDescription", {})
        name = ""
        for key in ("recommendedName", "submissionNames", "alternativeNames"):
            val = desc.get(key)
            if val:
                if isinstance(val, list):
                    name = val[0].get("fullName", {}).get("value", "")
                else:
                    name = val.get("fullName", {}).get("value", "")
                if name:
                    break

        results.append({
            "accession": entry["primaryAccession"],
            "name": name,
            "sequence": entry.get("sequence", {}).get("value", ""),
            "length": entry.get("sequence", {}).get("length", 0),
        })
    return results
```

---

## 2. IEDB MHC-II API -- Helper T-Cell Epitope Prediction

### 2.1 Endpoint

**URL:** `POST http://tools-cluster-interface.iedb.org/tools_api/mhcii/`
**Content-Type:** `application/x-www-form-urlencoded` (standard form POST)

### 2.2 Parameters

| Parameter | Required | Description |
|-----------|----------|-------------|
| `method` | Yes | Prediction method (see below) |
| `sequence_text` | Yes | Amino acid sequence (or FASTA with URI-encoded newlines) |
| `allele` | Yes | Comma-separated allele names |
| `length` | No | Peptide length(s), comma-separated. Default: 15. Use "asis" to predict as-is. |

### 2.3 Available Methods

| Method | Response Columns | Score Column | Notes |
|--------|-----------------|--------------|-------|
| `nn_align` | allele, seq_num, start, end, length, core_peptide, peptide, **ic50**, rank, adjusted_rank | `ic50` (nM) + `rank` (percentile) | Has `adjusted_rank` column |
| `NetMHCIIpan` | allele, seq_num, start, end, length, core_peptide, peptide, **ic50**, rank | `ic50` (nM) + `rank` (percentile) | Most allele coverage |
| `recommended` | allele, seq_num, start, end, length, core_peptide, peptide, **score**, rank | `score` (0-1, higher=better binder) + `rank` (percentile) | Consensus method |

**Key difference from MHC-I:** The MHC-II `recommended` method returns a `score`
column (0-1 scale, higher is better binding) instead of `ic50`. The `nn_align`
and `NetMHCIIpan` methods return `ic50` in nM (lower is better).

### 2.4 Canine DLA Alleles -- NOT SUPPORTED

**Tested and failed (all return "Invalid allele name"):**
- `DLA-DQA1*00101/DLA-DQB1*00201`
- `DLA-DQA1*001:01/DLA-DQB1*002:01`
- `DLA-DQA1*0101/DLA-DQB1*0201`
- `DLA-DQB1*00201`
- `DLA-DRB1*00101`
- `DLA-DRB1*01101`
- `DLA-DRB1*011:01`
- `DLA-DQA10101/DLA-DQB10201`

**Conclusion:** The IEDB MHC-II API does not support canine DLA class II alleles.

**Note on MHC-I vs MHC-II:** The existing pipeline (`04_immunogenicity.py`) uses
MHC-I with DLA-88 alleles (`DLA-8850101`, `DLA-8850801`, `DLA-8803401`), which
IS supported by IEDB MHC-I. The MHC-II API has a different allele catalog.

### 2.5 Fallback Strategy: HLA Alleles

Since canine DLA-II is not available, use HLA-DRB1 alleles as a cross-species proxy.
This is an accepted approach in veterinary vaccinology literature (DLA and HLA share
structural homology in the peptide-binding groove).

**Recommended HLA-II alleles for cross-species prediction:**

| Allele | Supertype | Rationale |
|--------|-----------|-----------|
| `HLA-DRB1*01:01` | DR1 | Most commonly used reference allele |
| `HLA-DRB1*04:01` | DR4 | Broad population coverage |
| `HLA-DRB1*07:01` | DR7 | Distinct binding specificity |
| `HLA-DRB1*11:01` | DR11 | Additional coverage |
| `HLA-DRB1*15:01` | DR15 | Additional coverage |

All tested successfully. Multiple alleles can be passed comma-separated in one request.

### 2.6 Response Format

Tab-separated values (TSV). Rows sorted by binding affinity (best binders first).

**Example (NetMHCIIpan, first 3 rows):**
```
allele	seq_num	start	end	length	core_peptide	peptide	ic50	rank
HLA-DRB1*01:01	1	34	48	15	INKLAGQPA	NITYLINKLAGQPAG	7.54	1.1
HLA-DRB1*01:01	1	35	49	15	INKLAGQPA	ITYLINKLAGQPAGD	8.71	1.6
```

**Binding thresholds (for ic50-based methods):**
- Strong binder: IC50 < 50 nM (rank < 2%)
- Moderate binder: IC50 < 500 nM (rank < 10%)
- Weak/non-binder: IC50 >= 500 nM

**Binding thresholds (for recommended method, score-based):**
- Strong binder: score > 0.7 (rank < 2%)
- Moderate binder: score > 0.3 (rank < 10%)

### 2.7 Timeout and Behavior

- Typical response time: 5-30 seconds depending on sequence length and allele count
- Use timeout of at least 120 seconds
- No rate-limit headers observed, but use polite delay (1-2s between calls)
- Multiple alleles in one call is supported and preferred over separate calls

### 2.8 Recommended Implementation

```python
import csv
import io
import time
from typing import Optional

import requests

IEDB_MHCII_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"

# HLA-DRB1 alleles as cross-species proxy for canine DLA class II
MHCII_ALLELES = [
    "HLA-DRB1*01:01",
    "HLA-DRB1*04:01",
    "HLA-DRB1*07:01",
]

MHCII_METHOD = "NetMHCIIpan"  # Best allele coverage; use "recommended" for consensus
MHCII_PEPTIDE_LENGTH = 15
MHCII_IC50_THRESHOLD = 500.0  # nM


def predict_mhcii_binding(
    sequence: str,
    alleles: Optional[list[str]] = None,
    method: str = MHCII_METHOD,
    length: int = MHCII_PEPTIDE_LENGTH,
) -> list[dict]:
    """Predict MHC-II binding using the IEDB API.

    Args:
        sequence: Amino acid sequence (minimum ~15 residues).
        alleles: List of HLA-DRB1 allele names. Defaults to MHCII_ALLELES.
        method: Prediction method. One of: NetMHCIIpan, nn_align, recommended.
        length: Peptide length for prediction (typically 15 for MHC-II).

    Returns:
        List of dicts with keys: allele, start, end, core_peptide, peptide,
        ic50 (or score for recommended method), rank.
    """
    if alleles is None:
        alleles = MHCII_ALLELES

    payload = {
        "method": method,
        "sequence_text": sequence,
        "allele": ",".join(alleles),
        "length": str(length),
    }

    resp = requests.post(IEDB_MHCII_URL, data=payload, timeout=120)
    resp.raise_for_status()

    return _parse_mhcii_response(resp.text, method)


def _parse_mhcii_response(text: str, method: str) -> list[dict]:
    """Parse TSV response from IEDB MHC-II API.

    Handles both ic50-based methods (nn_align, NetMHCIIpan) and
    score-based methods (recommended).
    """
    results = []

    # Check for error messages (no tab in first line)
    first_line = text.strip().split("\n")[0]
    if "\t" not in first_line:
        raise ValueError(f"IEDB API error: {first_line}")

    reader = csv.DictReader(io.StringIO(text), delimiter="\t")

    for row in reader:
        try:
            entry = {
                "allele": row["allele"].strip(),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "core_peptide": row["core_peptide"].strip(),
                "peptide": row["peptide"].strip(),
                "rank": float(row["rank"]),
            }

            # The score column differs by method
            if method == "recommended":
                entry["score"] = float(row["score"])
                entry["ic50"] = None
            else:
                entry["ic50"] = float(row["ic50"])
                entry["score"] = None

            results.append(entry)
        except (KeyError, ValueError):
            continue

    return results


def calculate_mhcii_immunogenicity(
    predictions: list[dict],
    method: str = MHCII_METHOD,
) -> float:
    """Score helper T-cell immunogenicity from MHC-II predictions.

    Returns fraction of peptides that are moderate-to-strong binders.
    """
    if not predictions:
        return 0.0

    best_per_peptide: dict[str, float] = {}

    for pred in predictions:
        peptide = pred["peptide"]

        if method == "recommended":
            # Higher score = better binder
            val = pred["score"]
            if peptide not in best_per_peptide or val > best_per_peptide[peptide]:
                best_per_peptide[peptide] = val
        else:
            # Lower IC50 = better binder
            val = pred["ic50"]
            if peptide not in best_per_peptide or val < best_per_peptide[peptide]:
                best_per_peptide[peptide] = val

    total = len(best_per_peptide)
    if total == 0:
        return 0.0

    if method == "recommended":
        good = sum(1 for s in best_per_peptide.values() if s > 0.3)
    else:
        good = sum(1 for ic50 in best_per_peptide.values() if ic50 < MHCII_IC50_THRESHOLD)

    return good / total
```

---

## 3. Integration Notes for the Marley Pipeline

### 3.1 Current State

The existing `pipeline/04_immunogenicity.py` uses:
- IEDB **MHC-I** API at `tools_api/mhci/` with `netmhcpan_ba` method
- DLA-88 class I alleles: `DLA-8850101`, `DLA-8850801`, `DLA-8803401`
- 9-mer peptides, IC50 threshold of 500 nM

### 3.2 What the MHC-II Integration Adds

MHC-II predicts **helper T-cell** (CD4+) epitopes, complementing MHC-I which
predicts **cytotoxic T-cell** (CD8+) epitopes. For vaccine design, both are needed:
- MHC-I: identifies peptides that trigger killer T cells
- MHC-II: identifies peptides that trigger helper T cells (critical for immune memory)

### 3.3 Recommended Architecture Changes

1. **Add UniProt fetcher module** -- `pipeline/utils/uniprot.py`
   - `fetch_sequence_by_accession(accession) -> str`
   - `search_sequences(name, organism_id) -> list[dict]`
   - Replace hardcoded FASTA files with dynamic fetching where needed

2. **Add MHC-II prediction** -- extend `pipeline/04_immunogenicity.py`
   - Add `predict_mhcii_binding()` alongside existing `predict_binding()` (MHC-I)
   - Compute separate `mhci_score` and `mhcii_score`
   - Combined immunogenicity = weighted average (e.g., 0.5 * mhci + 0.5 * mhcii)

3. **Accession mapping** -- create a curated mapping for target antigens:
   ```python
   ANTIGEN_ACCESSIONS = {
       "A2": "A4HZU7",           # A2 protein, L. infantum (found via search)
       "LACK": "J9XRQ3",         # Antigen LACK, L. infantum (found via search)
       "HSP70": "A0A6L0XI79",    # Hsp70, L. infantum, 654 aa (found via search)
       # KMP-11 and GP63 need manual curation -- not found by simple name search
       # LiHyp1 not in UniProt -- may need TriTrypDB accession
   }
   ```

### 3.4 Accession Warnings

The accessions originally proposed in the task do NOT correspond to the expected proteins:
- **Q25327** is a phosphatase from L. mexicana, not "A2 antigen"
- **P15711** is a Theileria parva protein, not "KMP-11"
- **Q25223** is an RNA polymerase fragment from L. amazonensis, not "GP63"
- **Q25292** is indeed an HSP70 fragment from L. infantum, but labeled as "LACK"

Use the search API or curate accessions from the literature before building a static mapping.
