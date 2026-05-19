#!/usr/bin/env python3
"""
DEPRECATED. This script targets the pre-Hydra `copal analyze --labels ...`
CLI and the old per-confidence `curated_conf{N}_no_leakage.csv` layout, both
of which have been superseded. It is retained for reference only and is not
maintained. For the current workflow see `docs/reproduce_paper.md`.

Extended COPAL analysis driver.

Like ``copal analyze`` (i.e. ``copal.analysis.pipeline.run_generate``), but with
broader output coverage:

  * Per-conf no-leakage variants for ``ranking_performance``, ``acyl_moiety``,
    ``chemical_similarity`` (each conf uses its own
    ``curated_conf{N}_no_leakage.csv``).
  * Per-model ``acyl_moiety`` splits emitted for conf 2 and conf 3 too.
  * ``mpkt`` scoped to RefAHL across conf {1, 2, 3} as one set
    (``mpkt_all_conf.csv``, ``mpkt_quartile_top6_all_conf.csv``,
    ``mpkt_correlation_all_conf.csv`` plus ``_no_leakage`` variants which use
    the union of the three curated sets).
  * Ablation tables (raw + summary, both regular and ``_no_leakage``) on
    RefAHL conf {1, 2, 3}, suffix-excluded; ``_no_leakage`` uses the union of
    the three curated sets.
  * ``full_rankings`` unchanged (all splits, no ``_no_leakage``).

The stock ``copal analyze`` CLI only emits ``_no_leakage`` for conf 1 and runs
ablation on conf 1 only.

Inputs read from ``--workflow-dir``:
  * ``per_model_rankings.csv``
  * ``ensemble_rankings.csv``
  * ``ablation_table.csv`` (loaded for parity but ablation is recomputed here)

Run from the pub repo root::

    .venv/bin/python applications/copal_analysis_extended.py \\
        --workflow-dir   <dir with per_model_rankings.csv etc.> \\
        --output-root    <where results/ should be written>
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from copal.analysis.pipeline import (
    AGGREGATIONS,
    DEFAULT_AGGREGATION,
    DEFAULT_RRF_K,
    MODEL_NAME_MAP,
    apply_global_suffix_exclusion,
    build_ablation_summary,
    build_ablation_table,
    build_acyl_moiety,
    build_chemical_similarity,
    build_full_rankings,
    build_mpkt,
    build_mpkt_correlation_test,
    build_mpkt_quartile_top6_success,
    build_ranked_excel,
    build_ranking_performance,
    build_true_rank_lookup,
    compute_avg_agreement_all,
    load_inputs,
    load_no_leakage_names,
    make_ensemble_ranking,
    normalize_inputs,
    write_results_csv,
)


REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
DEFAULT_PER_MODEL_SPLIT = ["AF3", "Boltz2", "Boltz2Affinity", "RF3"]
ABLATION_SUFFIX_EXCLUSIONS = ("SHORT", "MED", "LONG")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--workflow-dir", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument(
        "--labels-path",
        type=Path,
        default=DATA_DIR / "benchmark_data_all_corrected.csv",
    )
    parser.add_argument(
        "--ligands-path",
        type=Path,
        default=DATA_DIR / "Final_confirmed_ligands_58_annotated.csv",
    )
    parser.add_argument(
        "--no-leakage-conf1",
        type=Path,
        default=DATA_DIR / "curated_conf1_no_leakage.csv",
    )
    parser.add_argument(
        "--no-leakage-conf2",
        type=Path,
        default=DATA_DIR / "curated_conf2_no_leakage.csv",
    )
    parser.add_argument(
        "--no-leakage-conf3",
        type=Path,
        default=DATA_DIR / "curated_conf3_no_leakage.csv",
    )
    parser.add_argument(
        "--extra-model",
        action="append",
        default=[],
        metavar="NAME=DISPLAY",
        help="Register an extra model in MODEL_NAME_MAP (e.g. 'ralph=RALPH'). "
             "Repeatable. Use when running on a baseline whose model is not in "
             "the canonical AF3/Boltz2/Boltz2Affinity/RF3 set.",
    )
    parser.add_argument(
        "--per-model-split-models",
        type=str,
        default=",".join(DEFAULT_PER_MODEL_SPLIT),
        help="Comma-separated display names used for per-model acyl_moiety splits "
             f"(default: {','.join(DEFAULT_PER_MODEL_SPLIT)}).",
    )
    parser.add_argument(
        "--skip-ablation",
        action="store_true",
        help="Skip ablation table generation (use for single-model baselines like RALPH-BLAST).",
    )
    parser.add_argument(
        "--aggregation",
        type=str,
        default=DEFAULT_AGGREGATION,
        choices=AGGREGATIONS,
        help=(
            "Per-model rank aggregator used to re-derive the ensemble ranking "
            "from per_model_rankings.csv (overrides the ensemble_rankings.csv "
            "produced by `copal analyze`)."
        ),
    )
    parser.add_argument(
        "--rrf-k",
        type=int,
        default=DEFAULT_RRF_K,
        help=f"RRF damping constant (default {DEFAULT_RRF_K}; only used when --aggregation=rrf).",
    )
    return parser.parse_args()


def main() -> None:
    import sys as _sys

    print(
        "[DEPRECATED] applications/copal_analysis_extended.py is deprecated; "
        "use 'uv run copal analyze' (Hydra-backed) instead. "
        "See src/copal/config/set/*.yaml for available named sets.",
        file=_sys.stderr,
    )
    args = parse_args()
    workflow_dir = args.workflow_dir
    output_root = args.output_root
    ablation_path = workflow_dir / "ablation_table.csv"

    # Register any extra models (e.g. RALPH-BLAST baseline)
    from copal.analysis import pipeline as analysis_pipeline
    for spec in args.extra_model:
        if "=" not in spec:
            raise ValueError(f"--extra-model must be NAME=DISPLAY, got: {spec!r}")
        name, display = spec.split("=", 1)
        analysis_pipeline.MODEL_NAME_MAP[name.strip()] = display.strip()
        if display.strip() not in analysis_pipeline.RANKING_MODELS_ORDER:
            analysis_pipeline.RANKING_MODELS_ORDER = list(
                analysis_pipeline.RANKING_MODELS_ORDER
            ) + [display.strip()]
        print(f"Registered extra model: {name.strip()} -> {display.strip()}")

    per_model_split = [
        m.strip() for m in args.per_model_split_models.split(",") if m.strip()
    ]
    print(f"Per-model acyl_moiety split list: {per_model_split}")

    # If ablation table file is missing, fall back to a synthesized empty path
    # (run_generate's load_inputs requires a path; we recompute ablation below
    # anyway, so the loaded value is unused).
    if not ablation_path.exists():
        ablation_path = workflow_dir / "per_model_rankings.csv"  # any existing CSV

    labels_raw, rankings_raw, ensemble_raw, ligands_raw, _ = load_inputs(
        workflow_dir, args.ligands_path, ablation_path, args.labels_path,
    )
    labels, rankings, ensemble, ligands = normalize_inputs(
        labels_raw, rankings_raw, ensemble_raw, ligands_raw,
    )
    labels = apply_global_suffix_exclusion(labels)

    # Re-derive the ensemble from per_model_rankings under the requested
    # aggregation, regardless of what `copal analyze` originally wrote into
    # ensemble_rankings.csv. Lets callers sweep aggregators without re-running
    # postprocess. Schema matches the on-disk file.
    canonical_models = ["af3", "boltz2", "boltz2affinity", "rf3"]
    rerank_models = [m for m in canonical_models if m in rankings["model"].unique()]
    print(
        f"Re-deriving ensemble: aggregation={args.aggregation} "
        f"(rrf_k={args.rrf_k}) over {rerank_models}"
    )
    ensemble = make_ensemble_ranking(
        rankings, rerank_models, aggregation=args.aggregation, rrf_k=args.rrf_k,
    )
    # Mirror normalize_inputs: downstream consumers (build_chemical_similarity,
    # build_full_rankings, …) read `luxr` and `ahl` directly off the ensemble
    # frame, so split `luxr_ahl` here.
    _split = ensemble["luxr_ahl"].str.rsplit("_", n=1)
    ensemble["luxr"] = _split.str[0]
    ensemble["ahl"] = _split.str[1]

    nl_per_conf = {
        1: load_no_leakage_names(args.no_leakage_conf1),
        2: load_no_leakage_names(args.no_leakage_conf2),
        3: load_no_leakage_names(args.no_leakage_conf3),
    }
    nl_union = nl_per_conf[1] | nl_per_conf[2] | nl_per_conf[3]
    overlap = (
        len(nl_per_conf[1]) + len(nl_per_conf[2]) + len(nl_per_conf[3]) - len(nl_union)
    )
    print(
        f"no_leakage sizes: conf1={len(nl_per_conf[1])} "
        f"conf2={len(nl_per_conf[2])} conf3={len(nl_per_conf[3])} "
        f"union={len(nl_union)} overlap_lost={overlap}"
    )

    ranking_dir = output_root / "ranking_performance"
    acyl_dir = output_root / "acyl_moiety"
    mcs_dir = output_root / "chemical_similarity"
    mpkt_dir = output_root / "mpkt"
    full_dir = output_root / "full_rankings"
    ablation_dir = output_root / "ablation"
    for d in [ranking_dir, acyl_dir, mcs_dir, mpkt_dir, full_dir, ablation_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # ---- Per-conf outputs (regular + no_leakage) ----
    for conf in (1, 2, 3):
        base_subset = labels[
            (labels["split"].str.upper() == "REFAHL") & (labels["confidence"] == conf)
        ].drop_duplicates(subset=["LuxR_name", "split", "confidence"])
        print(f"[conf{conf}] base subset: {len(base_subset)} LuxRs")
        base_lookup = build_true_rank_lookup(base_subset, rankings, ensemble)

        write_results_csv(
            build_ranking_performance(base_subset, base_lookup),
            ranking_dir / f"ranking_performance_conf{conf}.csv",
        )
        write_results_csv(
            build_acyl_moiety(base_subset, ligands, base_lookup),
            acyl_dir / f"acyl_moiety_conf{conf}.csv",
        )
        for model_name in per_model_split:
            write_results_csv(
                build_acyl_moiety(
                    base_subset, ligands, base_lookup, rank_model_name=model_name
                ),
                acyl_dir / f"acyl_moiety_conf{conf}_{model_name}.csv",
            )
        write_results_csv(
            build_chemical_similarity(base_subset, ensemble, ligands),
            mcs_dir / f"chemical_similarity_conf{conf}.csv",
        )

        nl_subset = base_subset[
            base_subset["LuxR_name"].isin(nl_per_conf[conf])
        ].copy()
        print(f"[conf{conf}] no_leakage subset: {len(nl_subset)} LuxRs")
        if len(nl_subset) == 0:
            print(f"[conf{conf}] no_leakage subset empty -- skipping no_leakage outputs")
            continue
        nl_lookup = build_true_rank_lookup(nl_subset, rankings, ensemble)

        write_results_csv(
            build_ranking_performance(nl_subset, nl_lookup),
            ranking_dir / f"ranking_performance_conf{conf}_no_leakage.csv",
        )
        write_results_csv(
            build_acyl_moiety(nl_subset, ligands, nl_lookup),
            acyl_dir / f"acyl_moiety_conf{conf}_no_leakage.csv",
        )
        for model_name in per_model_split:
            write_results_csv(
                build_acyl_moiety(
                    nl_subset, ligands, nl_lookup, rank_model_name=model_name
                ),
                acyl_dir / f"acyl_moiety_conf{conf}_{model_name}_no_leakage.csv",
            )
        write_results_csv(
            build_chemical_similarity(nl_subset, ensemble, ligands),
            mcs_dir / f"chemical_similarity_conf{conf}_no_leakage.csv",
        )

    # ---- mpkt: RefAHL conf {1,2,3} as one set ----
    mpkt_subset = labels[
        (labels["split"].str.upper() == "REFAHL")
        & labels["confidence"].isin([1, 2, 3])
    ].drop_duplicates(subset=["LuxR_name", "split", "confidence"])
    print(f"[mpkt] base subset: {len(mpkt_subset)} LuxRs (RefAHL conf 1+2+3)")
    mpkt_lookup = build_true_rank_lookup(mpkt_subset, rankings, ensemble)
    ensemble_true = (
        mpkt_lookup[mpkt_lookup["model_name"] == "Ensemble"][["LuxR_name", "true_rank"]]
        .drop_duplicates()
    )
    mpkt_df = build_mpkt(mpkt_subset, rankings, ensemble_true)
    write_results_csv(mpkt_df, mpkt_dir / "mpkt_all_conf.csv")
    build_mpkt_quartile_top6_success(mpkt_df).to_csv(
        mpkt_dir / "mpkt_quartile_top6_all_conf.csv", index=False
    )
    build_mpkt_correlation_test(mpkt_df, mpkt_dir / "mpkt_correlation_all_conf.csv")

    nl_mpkt_subset = mpkt_subset[mpkt_subset["LuxR_name"].isin(nl_union)].copy()
    print(f"[mpkt] no_leakage(union) subset: {len(nl_mpkt_subset)} LuxRs")
    nl_mpkt_lookup = build_true_rank_lookup(nl_mpkt_subset, rankings, ensemble)
    nl_ensemble_true = (
        nl_mpkt_lookup[nl_mpkt_lookup["model_name"] == "Ensemble"][
            ["LuxR_name", "true_rank"]
        ]
        .drop_duplicates()
    )
    nl_mpkt_df = build_mpkt(nl_mpkt_subset, rankings, nl_ensemble_true)
    write_results_csv(nl_mpkt_df, mpkt_dir / "mpkt_all_conf_no_leakage.csv")
    build_mpkt_quartile_top6_success(nl_mpkt_df).to_csv(
        mpkt_dir / "mpkt_quartile_top6_all_conf_no_leakage.csv", index=False
    )
    build_mpkt_correlation_test(
        nl_mpkt_df, mpkt_dir / "mpkt_correlation_all_conf_no_leakage.csv"
    )

    if args.skip_ablation:
        print("--skip-ablation set; skipping ablation outputs")
    else:
        _run_ablation(
            labels, rankings, nl_union, ablation_dir,
            aggregation=args.aggregation, rrf_k=args.rrf_k,
        )

    # ---- full_rankings: unchanged (all splits) ----
    full_rankings_df = build_full_rankings(labels, rankings, ensemble, ligands)
    write_results_csv(full_rankings_df, full_dir / "all_splits.csv")

    rank_pivot_all = (
        rankings.pivot_table(
            index=["luxr", "ahl"], columns="model", values="rank", aggfunc="first"
        )
        .reset_index()
    )
    model_rank_cols_all = [m for m in MODEL_NAME_MAP if m in rank_pivot_all.columns]
    mpkt_all_recs = []
    for luxr, group in rank_pivot_all.groupby("luxr"):
        mpkt_all_recs.append(
            {
                "LuxR_name": luxr,
                "mpKT": compute_avg_agreement_all(group, model_rank_cols_all),
            }
        )
    mpkt_all_df = pd.DataFrame(mpkt_all_recs)

    excel_models = [
        ("AF3", "af3"),
        ("Boltz2", "boltz2"),
        ("Boltz2Affinity", "boltz2affinity"),
        ("RF3", "rf3"),
        ("ENSEMBLE", None),
    ]
    # Add any extra models registered via --extra-model
    seen = {key for _, key in excel_models}
    for spec in args.extra_model:
        name, display = spec.split("=", 1)
        if name.strip() not in seen:
            excel_models.append((display.strip(), name.strip()))

    for display_name, model_key in excel_models:
        excel_path = full_dir / f"{display_name}.xlsx"
        # Skip a model if its key isn't actually in rankings (avoids empty workbooks)
        if model_key is not None and model_key not in rankings["model"].unique():
            continue
        build_ranked_excel(
            model_key, rankings, ensemble, labels, ligands, mpkt_all_df, excel_path
        )
        print(f"Wrote {excel_path}")

    print("DONE")


def _sort_ablation_by_inclusion(table):
    # Top row: all four models. Bottom row: af3 only.
    # Primary: model count descending. Secondary: [af3, boltz2, boltz2 affinity, rf3] ascending,
    # so within each count group False sorts before True on each flag — placing the af3-only
    # row last among singletons.
    if table.empty:
        return table
    out = table.copy()
    out["_n"] = (
        out["af3 included"].astype(int)
        + out["boltz2 included"].astype(int)
        + out["boltz2 affinity included"].astype(int)
        + out["rf3 included"].astype(int)
    )
    out = out.sort_values(
        ["_n", "af3 included", "boltz2 included", "boltz2 affinity included", "rf3 included"],
        ascending=[False, True, True, True, True],
    ).reset_index(drop=True)
    return out.drop(columns="_n")


def _run_ablation(
    labels, rankings, nl_union, ablation_dir,
    aggregation=DEFAULT_AGGREGATION, rrf_k=DEFAULT_RRF_K,
):
    # ---- Ablation: RefAHL conf {1,2,3}, suffix-excluded ----
    abl_pool = labels[
        labels["true_ahl_id"].notna()
        & (labels["split"].str.upper() == "REFAHL")
        & labels["confidence"].isin([1, 2, 3])
    ].copy()
    abl_pool = abl_pool[
        ~abl_pool["LuxR_name"].str.endswith(ABLATION_SUFFIX_EXCLUSIONS, na=False)
    ]
    abl_pool["luxr"] = abl_pool["LuxR_name"]
    print(f"[ablation] full pool: {len(abl_pool)} LuxRs")

    abl_table = build_ablation_table(
        rankings, abl_pool, aggregation=aggregation, rrf_k=rrf_k,
    )
    _sort_ablation_by_inclusion(abl_table).to_csv(
        ablation_dir / "ablation_table_raw.csv", index=False
    )
    build_ablation_summary(abl_table).to_csv(
        ablation_dir / "ablation_table.csv", index=False
    )

    nl_abl_pool = abl_pool[abl_pool["LuxR_name"].isin(nl_union)].copy()
    print(f"[ablation] no_leakage pool: {len(nl_abl_pool)} LuxRs")
    nl_abl_table = build_ablation_table(
        rankings, nl_abl_pool, aggregation=aggregation, rrf_k=rrf_k,
    )
    _sort_ablation_by_inclusion(nl_abl_table).to_csv(
        ablation_dir / "ablation_table_raw_no_leakage.csv", index=False
    )
    build_ablation_summary(nl_abl_table).to_csv(
        ablation_dir / "ablation_table_no_leakage.csv", index=False
    )

    sorted_cols = [
        "top 5 success rate", "top 6 success rate", "top 3 success rate",
        "top 1 success rate", "af3 included", "boltz2 included",
        "boltz2 affinity included", "rf3 included",
    ]
    nl_sorted = nl_abl_table.copy()
    nl_sorted["_n"] = (
        nl_sorted["af3 included"].astype(int)
        + nl_sorted["boltz2 included"].astype(int)
        + nl_sorted["boltz2 affinity included"].astype(int)
        + nl_sorted["rf3 included"].astype(int)
    )
    nl_sorted = nl_sorted.sort_values(
        ["_n", "af3 included", "boltz2 included", "boltz2 affinity included", "rf3 included"],
        ascending=[False, False, False, False, False],
    ).reset_index(drop=True)
    nl_sorted[sorted_cols].to_csv(
        ablation_dir / "sorted_ablation_full_no_leakage.csv", index=False
    )


if __name__ == "__main__":
    main()
