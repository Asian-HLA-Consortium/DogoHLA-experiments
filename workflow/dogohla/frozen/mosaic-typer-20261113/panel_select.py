"""Panel selection for the proposer with a hard leave-donor-and-family-out check.

Arm C: P = IPD genomic alleles of the gene only.
Arm D: P = IPD + the cohort cores that the arm-B reference of the donor's fold carries (its PG
       candidates); every cohort row must be carried by at least one donor outside the donor's
       FOLDS.json component (component from atlas panels.build_exclusion_map_from_folds).
"""
from __future__ import annotations

from typing import Iterable, List, Optional, Sequence, Set


class LeakageError(AssertionError):
    pass


def check_no_leakage(rows: Iterable[int], row_is_ipd: Sequence[bool],
                     row_carriers: Sequence[Set[str]], excluded: Set[str], donor: str,
                     component_absent_rows: Set[int], arm: str) -> None:
    """Raise LeakageError if any panel row is justified only by the donor or its relatives.

    row_is_ipd[r]        row r is a public IPD genomic allele
    row_carriers[r]      cohort donors whose truth core equals row r
    component_absent_rows rows equal to an IPD-absent truth core of a donor in the component
    """
    if donor not in excluded:
        raise LeakageError('donor missing from its own component')
    rows = list(rows)
    for r in rows:
        if row_is_ipd[r]:
            continue
        if arm == 'C':
            raise LeakageError(f'arm C panel holds a non-IPD row {r}')
        if not (row_carriers[r] - excluded):
            raise LeakageError(f'row {r} is carried only by the donor/component {sorted(row_carriers[r])}')
    if arm == 'C':
        bad = component_absent_rows & set(rows)
        if bad:
            raise LeakageError(f'IPD-absent truth of the component present in arm C panel: {sorted(bad)}')
    else:
        for r in component_absent_rows & set(rows):
            if not (row_carriers[r] - excluded):
                raise LeakageError(f'component truth row {r} has no unrelated carrier')


def select_rows(arm: str, row_is_ipd: Sequence[bool], row_seq_sha: Sequence[str],
                pg_shas: Optional[Set[str]] = None) -> List[int]:
    if arm == 'C':
        return [r for r, ipd in enumerate(row_is_ipd) if ipd]
    if arm == 'D':
        assert pg_shas is not None
        return [r for r, ipd in enumerate(row_is_ipd) if ipd or row_seq_sha[r] in pg_shas]
    raise ValueError(arm)
