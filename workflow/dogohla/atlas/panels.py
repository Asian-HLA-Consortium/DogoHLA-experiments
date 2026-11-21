"""Leave-donor-and-family-out panel construction and input loaders.

Uses $S/inputs/FOLDS.json (union-find over its 'groups' components and 'donors' registry)
as well as simulation-augmented-regional-20261024/fold0..fold4/references_core/COMPLETE.json.
Enforces hard rule 2: the held-out donor and its relatives never contribute to P(q).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

GENES: Tuple[str, ...] = ("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
CLASS_I_GENES: Set[str] = {"A", "B", "C"}
CLASS_II_GENES: Set[str] = {"DPA1", "DPB1", "DQA1", "DQB1", "DRB1"}
FLANK: int = 2000


def sha256_str(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def load_source_records(path: Path) -> List[Dict]:
    data = json.loads(Path(path).read_text())
    records = data["records"]
    for r in records:
        core = r["sequence"][FLANK:-FLANK]
        if sha256_str(core) != r["core_sha256"]:
            raise ValueError(f"Core sha256 mismatch for {r.get('source_id')}")
        r["core_sequence"] = core
        r["core_exons"] = [(s - FLANK, e - FLANK) for s, e in r["exons"]]
    return records


def load_ipd_reference(fasta_path: Path) -> Dict[str, List[Dict]]:
    """Load IPD genomic alleles per classical gene from genomic.fa.

    Header format in genomic.fa:
      >HLA-A*01:01:01:01 <num_exons> <start_0> <end_0_inclusive> ...
    """
    by_gene: Dict[str, List[Dict]] = {g: [] for g in GENES}
    name: Optional[str] = None
    header_ints: List[int] = []
    seq_parts: List[str] = []

    def _flush(allele_name: Optional[str], ints: List[int], parts: List[str]) -> None:
        if not allele_name:
            return
        gene = allele_name.split("*")[0].removeprefix("HLA-")
        if gene not in by_gene:
            return
        seq = "".join(parts).upper()
        n_exons = ints[0] if ints else 0
        exons = [(ints[1 + 2 * i], ints[2 + 2 * i] + 1) for i in range(n_exons)]
        by_gene[gene].append(
            {
                "allele": allele_name,
                "gene": gene,
                "sequence": seq,
                "exons": exons,
                "sequence_sha256": sha256_str(seq),
            }
        )

    with Path(fasta_path).open() as f:
        for line in f:
            if line.startswith(">"):
                _flush(name, header_ints, seq_parts)
                tokens = line[1:].strip().split()
                name = tokens[0]
                header_ints = [int(x) for x in tokens[1:]]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        _flush(name, header_ints, seq_parts)
    return by_gene


def build_exclusion_map_from_folds(
    folds_source: Union[Path, str, Dict],
    truth_donors: Optional[Iterable[str]] = None,
) -> Dict[str, Set[str]]:
    """Build donor -> excluded_donors map via union-find over FOLDS.json 'groups' components.

    Every donor in FOLDS.json ('groups' and 'donors') is placed in its connected component
    under 'groups'. Asserts that every donor in truth_donors is resolved in FOLDS.json.
    """
    if isinstance(folds_source, (Path, str)):
        folds_data = json.loads(Path(folds_source).read_text())
    else:
        folds_data = folds_source

    if "groups" not in folds_data:
        raise ValueError("FOLDS.json is missing required 'groups' field")

    groups: Sequence[Sequence[str]] = folds_data["groups"]
    registered_donors: Set[str] = set()
    if "donors" in folds_data:
        if isinstance(folds_data["donors"], dict):
            registered_donors.update(folds_data["donors"].keys())
        else:
            registered_donors.update(folds_data["donors"])
    for g in groups:
        for d in g:
            registered_donors.add(d)

    if truth_donors is not None:
        truth_set = set(truth_donors)
        unresolved = sorted(truth_set - registered_donors)
        if unresolved:
            raise AssertionError(
                f"Truth donors not resolved in FOLDS.json: {unresolved}"
            )

    parent: Dict[str, str] = {d: d for d in registered_donors}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for g in groups:
        members = [d for d in g if d in parent]
        if len(members) > 1:
            root0 = find(members[0])
            for m in members[1:]:
                parent[find(m)] = root0

    components: Dict[str, Set[str]] = {}
    for d in registered_donors:
        r = find(d)
        components.setdefault(r, set()).add(d)

    excl_map = {d: set(components[find(d)]) for d in registered_donors}
    if truth_donors is not None:
        for d in truth_donors:
            assert d in excl_map and d in excl_map[d], f"Unresolved donor {d}"
    return excl_map


def load_augmented_fold_grouping(
    augmented_root: Path,
) -> Tuple[Dict[str, int], Dict[int, Set[str]], Dict[int, Set[str]]]:
    """Load the 5-fold CV grouping from simulation-augmented-regional-20261024/fold0..fold4."""
    donor_to_fold: Dict[str, int] = {}
    fold_excluded: Dict[int, Set[str]] = {}
    fold_training: Dict[int, Set[str]] = {}

    for k in range(5):
        comp_path = Path(augmented_root) / f"fold{k}" / "references_core" / "COMPLETE.json"
        comp = json.loads(comp_path.read_text())
        if comp.get("status") != "complete":
            raise ValueError(f"Incomplete fold manifest: {comp_path}")
        excl = set(comp["excluded_test_donors"])
        train = set(comp["training_donors"])
        if excl & train:
            raise ValueError(f"Overlap between train and test in fold {k}")
        fold_excluded[k] = excl
        fold_training[k] = train
        for d in excl:
            if d in donor_to_fold:
                raise ValueError(f"Donor {d} in multiple test folds")
            donor_to_fold[d] = k

    return donor_to_fold, fold_excluded, fold_training


def verify_no_leakage(
    query_donor: str,
    excluded_donors: Set[str],
    panel_records: Iterable[Dict],
) -> None:
    """Hard rule 2: fail immediately if query_donor or any relative is in the panel."""
    if query_donor not in excluded_donors:
        raise ValueError(f"Query donor {query_donor} missing from its own excluded_donors set")
    leaked = sorted(
        {
            r["donor"]
            for r in panel_records
            if "donor" in r and r["donor"] is not None and r["donor"] in excluded_donors
        }
    )
    if leaked:
        raise ValueError(
            f"Held-out donor or family member leaked into panel for {query_donor}: {leaked}"
        )


def select_cohort_panel_for_donor(
    query_donor: str,
    gene_cohort_records: List[Dict],
    excluded_donors_map: Dict[str, Set[str]],
    allowed_resources: Optional[Set[str]] = None,
) -> List[Dict]:
    """Return cohort records outside query_donor's leave-donor-and-family-out group."""
    if query_donor not in excluded_donors_map:
        raise KeyError(f"Unknown donor {query_donor} in excluded_donors_map")
    excluded = excluded_donors_map[query_donor]
    panel = [
        r
        for r in gene_cohort_records
        if r.get("eligible_reference", True)
        and r["donor"] not in excluded
        and (allowed_resources is None or r["resource"] in allowed_resources)
    ]
    verify_no_leakage(query_donor, excluded, panel)
    return panel
