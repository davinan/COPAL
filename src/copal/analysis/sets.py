"""Benchmark set abstraction for copal analysis.

A `BenchmarkSet` is a resolved view of the labels CSV — already filtered by
split, confidence tier, no-leakage curation, and suffix exclusion. Downstream
analysis (`pipeline.run_analysis_for_set`) operates on a single resolved set,
without any internal looping over conf{1,2,3}/no_leakage variants.

Every filter axis is independently optional. An empty value (``[]``, ``None``,
or a missing YAML key) means "do not filter on that axis" — not "filter to
empty".
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from copal.analysis import pipeline as _pipeline


@dataclass(frozen=True)
class BenchmarkSetSpec:
    """Unresolved description of which rows belong to a benchmark set."""

    name: str
    splits: tuple[str, ...] = ()
    confidences: tuple[int, ...] = ()
    no_leakage_csvs: tuple[Path, ...] = ()
    exclude_suffixes: tuple[str, ...] = ()
    include_names: tuple[str, ...] = ()


@dataclass(frozen=True)
class BenchmarkSet:
    """A resolved benchmark set: post-filter labels + ligands metadata."""

    name: str
    labels: pd.DataFrame
    ligands: pd.DataFrame
    name_lookup: frozenset[str]
    provenance: dict[str, Any] = field(default_factory=dict)

    def is_empty(self) -> bool:
        return self.labels.empty

    def n_luxrs(self) -> int:
        if self.labels.empty:
            return 0
        return int(self.labels["LuxR_name"].nunique())


def load_no_leakage_union(paths: Iterable[Path]) -> frozenset[str]:
    """Union of `refahl_id` values across one or more curated no-leakage CSVs."""
    union: set[str] = set()
    for p in paths:
        union |= _pipeline.load_no_leakage_names(Path(p))
    return frozenset(union)


def _apply_split_filter(
    labels: pd.DataFrame, splits: Iterable[str]
) -> pd.DataFrame:
    targets = [s for s in (splits or ()) if s]
    if not targets:
        return labels
    upper = labels["split"].astype(str).str.upper()
    keep = {s.upper() for s in targets}
    return labels[upper.isin(keep)]


def _apply_confidence_filter(
    labels: pd.DataFrame, confs: Iterable[int]
) -> pd.DataFrame:
    targets = list(confs or ())
    if not targets:
        return labels
    return labels[labels["confidence"].isin(targets)]


def _apply_leakage_filter(
    labels: pd.DataFrame, sources: Iterable[Path]
) -> pd.DataFrame:
    paths = [Path(p) for p in (sources or ())]
    if not paths:
        return labels
    keep = load_no_leakage_union(paths)
    return labels[labels["LuxR_name"].isin(keep)]


def _apply_include_names_filter(
    labels: pd.DataFrame, names: Iterable[str]
) -> pd.DataFrame:
    targets = [n for n in (names or ()) if n]
    if not targets:
        return labels
    keep = {str(n).upper() for n in targets}
    upper = labels["LuxR_name"].astype(str).str.upper()
    return labels[upper.isin(keep)]


def _apply_suffix_filter(
    labels: pd.DataFrame, suffixes: Iterable[str]
) -> pd.DataFrame:
    targets = [s for s in (suffixes or ()) if s]
    if not targets:
        return labels
    upper = labels["LuxR_name"].astype(str).str.upper()
    mask = pd.Series(False, index=labels.index)
    for s in targets:
        mask |= upper.str.endswith(s.upper())
    return labels[~mask]


def _normalize_labels_only(labels_raw: pd.DataFrame) -> pd.DataFrame:
    """Apply the labels-half of `pipeline.normalize_inputs` without touching rankings/ligands."""
    labels_norm = labels_raw.copy()
    labels_norm["LuxR_name"] = labels_norm["refahl_id"].apply(_pipeline.normalize_token)
    labels_norm["split"] = labels_norm["split"].astype(str)
    labels_norm["confidence"] = pd.to_numeric(
        labels_norm["refahl_confidence_rating"], errors="coerce"
    )
    labels_norm["true_ahl_id"] = labels_norm["true_ahl_id"].apply(
        _pipeline.normalize_optional_token
    )
    return labels_norm[
        ["LuxR_name", "split", "confidence", "true_ahl_id", "true_ahl_abbrev"]
    ].drop_duplicates()


def _normalize_ligands(ligands_raw: pd.DataFrame) -> pd.DataFrame:
    """Apply the ligands-half of `pipeline.normalize_inputs`."""
    ligands_norm = ligands_raw.copy()
    ligands_norm["AHLid"] = ligands_norm["AHLid"].apply(_pipeline.normalize_token)
    return ligands_norm


def load_set(
    spec: BenchmarkSetSpec,
    labels_path: Path | None,
    ligands_path: Path | None,
) -> BenchmarkSet:
    """Resolve a `BenchmarkSetSpec` into a `BenchmarkSet` by reading + filtering CSVs.

    Either path may be ``None`` to signal inference mode for that axis: the
    corresponding DataFrame is returned empty with the expected columns, and
    the caller (typically `run_full_ranking_for_set`) synthesizes a minimal
    view from the rankings themselves.
    """
    if labels_path is None:
        labels = pd.DataFrame(
            columns=["LuxR_name", "split", "confidence", "true_ahl_id", "true_ahl_abbrev"]
        )
    else:
        labels = _normalize_labels_only(pd.read_csv(labels_path))
        labels = _apply_suffix_filter(labels, spec.exclude_suffixes)
        labels = _apply_split_filter(labels, spec.splits)
        labels = _apply_confidence_filter(labels, spec.confidences)
        labels = _apply_leakage_filter(labels, spec.no_leakage_csvs)
        labels = _apply_include_names_filter(labels, spec.include_names)
        labels = labels.drop_duplicates(
            subset=["LuxR_name", "split", "confidence"]
        ).reset_index(drop=True)

    if ligands_path is None:
        ligands = pd.DataFrame(columns=["AHLid", "abbreviation", "subcategory"])
    else:
        ligands = _normalize_ligands(pd.read_csv(ligands_path))

    name_lookup = frozenset(labels["LuxR_name"].astype(str).tolist())
    provenance: dict[str, Any] = {
        "splits": list(spec.splits),
        "confidences": list(spec.confidences),
        "no_leakage_csvs": [str(p) for p in spec.no_leakage_csvs],
        "exclude_suffixes": list(spec.exclude_suffixes),
        "include_names": list(spec.include_names),
        "labels_path": str(labels_path) if labels_path is not None else None,
        "ligands_path": str(ligands_path) if ligands_path is not None else None,
        "n_rows": int(len(labels)),
        "n_luxrs": int(labels["LuxR_name"].nunique()) if not labels.empty else 0,
    }
    return BenchmarkSet(
        name=spec.name,
        labels=labels,
        ligands=ligands,
        name_lookup=name_lookup,
        provenance=provenance,
    )


def from_omega(cfg_set) -> BenchmarkSetSpec:
    """Translate an OmegaConf `set` group into a `BenchmarkSetSpec`.

    Treats missing keys, ``[]``, and ``null`` interchangeably as "empty" (no
    filter on that axis).
    """
    raw_filter = cfg_set.get("filter") if hasattr(cfg_set, "get") else None
    if raw_filter is None:
        raw_filter = {}

    def _seq(key: str) -> tuple:
        v = raw_filter.get(key) if hasattr(raw_filter, "get") else None
        if not v:
            return ()
        return tuple(v)

    no_leakage_raw = raw_filter.get("no_leakage_csvs") if hasattr(raw_filter, "get") else None
    no_leakage = tuple(Path(p) for p in (no_leakage_raw or ()))

    confidences_raw = raw_filter.get("confidences") if hasattr(raw_filter, "get") else None
    confidences = tuple(int(c) for c in (confidences_raw or ()))

    return BenchmarkSetSpec(
        name=str(cfg_set.name),
        splits=_seq("splits"),
        confidences=confidences,
        no_leakage_csvs=no_leakage,
        exclude_suffixes=_seq("exclude_suffixes"),
        include_names=_seq("include_names"),
    )
