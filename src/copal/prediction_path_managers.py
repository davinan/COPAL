"""Per-model path managers used to discover, validate, and stage prediction
artifacts (structures + confidence files) produced by AF3, Boltz2, and RF3.
"""
from __future__ import annotations

import json
from abc import ABC, abstractmethod
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np
import pandas as pd


class ModelPathManager(ABC):
    def __init__(self, basepath: Path):
        self.basepath = Path(basepath)

    @abstractmethod
    def get_structure_files(self, chunk: str, name: str): ...

    @abstractmethod
    def get_confidence_files(self, chunk: str, name: str): ...

    @abstractmethod
    def get_deletable(self, chunk: str, name: str): ...

    def get(self, chunk, name):
        sfiles = self.get_structure_files(chunk, name)
        try:
            cfiles = self.get_confidence_files(chunk, name)
        except Exception:
            return False, {}

        deletable = self.get_deletable(chunk, name)

        sum_size_mb = lambda xs: sum(x.stat().st_size / (1024 * 1024) for x in xs)

        output = {
            "structures_to_save": sfiles["path_list"],
            "confidence_to_save": cfiles["path_list"],
            "to_delete": deletable,
            "structures_size_total_MB": sum_size_mb(sfiles["path_list"]),
            "confidence_size_total_MB": sum_size_mb(cfiles["path_list"]),
            "to_delete_total_MB": sum_size_mb(deletable),
        }
        output["saving_total_MB"] = (
            output["structures_size_total_MB"] + output["confidence_size_total_MB"]
        )
        return sfiles["has_all"] and cfiles["has_all"], output


class Boltz2PathManager(ModelPathManager):
    def pred_dir_iter(self):
        for chunk in self.basepath.iterdir():
            preds = chunk / f"boltz_results_{chunk.name}" / "predictions"
            if not preds.exists():
                continue
            for pred in preds.iterdir():
                yield chunk.name, pred.name

    def get_pred_dir(self, chunk, name):
        return self.basepath / str(chunk) / f"boltz_results_{chunk}" / "predictions" / name

    def get_structure_files(self, chunk, name):
        pdir = self.get_pred_dir(chunk, name)
        files = list(pdir.glob("*.cif"))
        return {"path_list": files, "has_all": len(files) == 5}

    def get_confidence_files(self, chunk, name):
        sample_dir = self.get_pred_dir(chunk, name)
        chain_i, chain_j = "A", "B"
        recs = []
        try:
            aa = strucio.load_structure(list(sample_dir.glob("*.cif"))[0])
        except Exception:
            return {"has_all": False, "path_list": []}

        chain_ids = []
        for residue in struc.residue_iter(aa):
            if residue.hetero[0]:
                chain_ids.extend(list(residue.chain_id))
            else:
                chain_ids.append(residue.chain_id[0])
        chain_ids = np.array(chain_ids)

        for pae_path in sample_dir.glob("pae_*.npz"):
            _, *name_split, __, ___ = pae_path.name.split("_")
            split_name = "_".join(name_split)
            model_index = pae_path.name.split("_")[-1].split(".npz")[0]

            is_interface = chain_ids[:, None] != chain_ids[None, :]
            is_chain_of_interest = (chain_ids == chain_i) | (chain_ids == chain_j)
            is_interface_of_interest = is_chain_of_interest[:, None] | is_chain_of_interest[None, :]
            mask = is_interface & is_interface_of_interest

            pae = np.load(pae_path)["pae"]
            conf_path = pae_path.parent / f"confidence_{split_name}_model_{model_index}.json"
            with open(conf_path) as fh:
                data = json.load(fh)

            recs.append({
                "name": split_name,
                "model": "boltz2",
                "batch_idx": model_index,
                "min_pae_interface": pae[mask].min(),
                "pae_interface": pae[mask].mean(),
                "iptm": data["iptm"],
            })

        for affinity_path in sample_dir.glob("affinity_*.json"):
            _, *name_split = affinity_path.name.split("_")
            name_split[-1] = name_split[-1].split(".json")[0]
            split_name = "_".join(name_split)
            with open(affinity_path) as fh:
                data = json.load(fh)
            recs.append({
                "name": split_name,
                "model": "boltz2affinity",
                "batch_idx": 0,
                "affinity_pred_value": data["affinity_pred_value"],
                "affinity_pred_value1": data["affinity_pred_value1"],
                "affinity_pred_value2": data["affinity_pred_value2"],
                "affinity_probability_binary": data["affinity_probability_binary"],
                "affinity_probability_binary2": data["affinity_probability_binary2"],
            })

        df = pd.DataFrame(recs)
        df.to_csv(sample_dir / "confidences.csv")
        files = list(sample_dir.glob("*.json")) + list(sample_dir.glob("*.csv"))
        return {"path_list": files, "has_all": len(files) >= 6}

    def get_deletable(self, chunk, name):
        return list(self.get_pred_dir(chunk, name).glob("*.npz"))


class AF3PathManager(ModelPathManager):
    def get_pred_dir(self, chunk, name):
        return self.basepath / str(chunk) / name

    def get_structure_files(self, chunk, name):
        pdir = self.get_pred_dir(chunk, name)
        files = list(pdir.glob("*/*.cif"))
        return {"path_list": files, "has_all": len(files) == 5}

    def get_confidence_files(self, chunk, name):
        base_dir = self.get_pred_dir(chunk, name)
        chain_i, chain_j = "A", "B"
        recs = []
        for sample_dir in base_dir.iterdir():
            if not sample_dir.is_dir() or "sample" not in sample_dir.name:
                continue
            batch_idx = int(sample_dir.name.split("-")[-1])
            conf_paths = list(sample_dir.glob("*confidences.json"))
            conf_path = [x for x in conf_paths if "summary" not in x.name]
            summ_path = [x for x in conf_paths if "summary" in x.name]
            if len(conf_path) != 1 or len(summ_path) != 1:
                continue

            with open(conf_path[0]) as fh:
                af3_conf = json.load(fh)
            pae = np.array(af3_conf["pae"])
            chain_id = np.array(af3_conf["token_chain_ids"])
            is_chain_of_interest = (chain_id == chain_i) | (chain_id == chain_j)
            is_interface_of_interest = is_chain_of_interest[:, None] | is_chain_of_interest[None, :]
            is_interface = chain_id[:, None] != chain_id[None, :]
            mask = is_interface & is_interface_of_interest

            with open(summ_path[0]) as fh:
                data = json.load(fh)

            recs.append({
                "seed": sample_dir.name.split("_")[0],
                "model": "af3",
                "batch_idx": batch_idx,
                "min_pae_interface": pae[mask].min(),
                "pae_interface": pae[mask].mean(),
                "iptm": data["chain_pair_iptm"][0][1],
                "min_pae_interface_af3": data["chain_pair_pae_min"][0][1],
            })

        agg_conf = base_dir / "confidences.csv"
        pd.DataFrame(recs).to_csv(agg_conf)

        files = list(base_dir.glob("*/*_summary_confidences.json")) + [agg_conf]
        files = [x for x in files if x.exists()]
        return {"path_list": files, "has_all": len(files) == 6}

    def get_deletable(self, chunk, name):
        base_dir = self.get_pred_dir(chunk, name)
        deletable = list(base_dir.glob("*.cif"))
        for sample_dir in base_dir.iterdir():
            if sample_dir.is_dir() and "sample" in sample_dir.name:
                deletable.append(sample_dir / f"{name}_{sample_dir.name}_confidences.json")
        return [x for x in deletable if x.exists()]

    def pred_dir_iter(self):
        for chunk in self.basepath.iterdir():
            if not chunk.is_dir():
                continue
            for pred in chunk.iterdir():
                if pred.is_dir():
                    yield chunk.name, pred.name


class RF3PathManager(ModelPathManager):
    """Handles the current production-branch RF3 output layout:

        <chunk>/<name>/
            <name>_ranking_scores.csv           (aggregate ranking per sample)
            <name>_summary_confidences.json     (aggregate summary)
            seed-*_sample-*/
                <name>_seed-*_sample-*_model.cif
                <name>_seed-*_sample-*_summary_confidences.json
                <name>_seed-*_sample-*_confidences.json

    Older RF3 output shape (top-level `.cif.gz` + a single `.score` CSV +
    `_metrics.csv`) is synthesized from the new shape on first touch, so
    `copal analyze` — which was written for the old shape — can keep
    reading the same `.score`/`_metrics.csv` schema unchanged.
    """

    def pred_dir_iter(self):
        for chunk in self.basepath.iterdir():
            if not chunk.is_dir():
                continue
            for pred in chunk.iterdir():
                if pred.is_dir() and any(pred.glob("seed-*_sample-*")):
                    yield chunk.name, pred.name

    def get_pred_dir(self, chunk, name):
        return self.basepath / str(chunk) / name

    def get_structure_files(self, chunk, name):
        pdir = self.get_pred_dir(chunk, name)
        files = sorted(pdir.glob("seed-*_sample-*/*_model.cif"))
        return {"path_list": files, "has_all": len(files) == 5}

    def _write_compat_csvs(self, pdir: Path, name: str) -> None:
        """Derive <name>.score + <name>_metrics.csv from per-sample summary JSONs."""
        score_path = pdir / f"{name}.score"
        metrics_path = pdir / f"{name}_metrics.csv"

        rows = []
        iptm_by_batch: dict[int, float] = {}
        for sample_dir in sorted(pdir.glob("seed-*_sample-*")):
            batch_idx = int(sample_dir.name.split("sample-")[-1])
            summ_path = sample_dir / f"{name}_{sample_dir.name}_summary_confidences.json"
            if not summ_path.exists():
                continue
            with open(summ_path) as fh:
                data = json.load(fh)
            iptm = data.get("iptm")
            # chain_pair_pae[0][1] is the A-B mean PAE; *_min is the A-B min PAE.
            pae_ab = _safe_pair_value(data.get("chain_pair_pae"))
            min_pae_ab = _safe_pair_value(data.get("chain_pair_pae_min"))
            rows.append({
                "batch_idx": batch_idx,
                "chain_i_interface": "A_1",
                "chain_j_interface": "B_1",
                "iptm": iptm,
                "pae_interface": pae_ab,
                "min_pae_interface": min_pae_ab,
            })
            if iptm is not None:
                iptm_by_batch[batch_idx] = float(iptm)

        if rows:
            pd.DataFrame(rows).to_csv(score_path, index=False)
            metrics_row = {f"iptm.iptm_protein_ligand_{k}": iptm_by_batch.get(k) for k in range(5)}
            pd.DataFrame([metrics_row]).to_csv(metrics_path, index=False)

    def get_confidence_files(self, chunk, name):
        pdir = self.get_pred_dir(chunk, name)
        score_path = pdir / f"{name}.score"
        metrics_path = pdir / f"{name}_metrics.csv"
        if not (score_path.exists() and metrics_path.exists()):
            self._write_compat_csvs(pdir, name)
        summary_jsons = list(pdir.glob("seed-*_sample-*/*_summary_confidences.json"))
        files = [p for p in [score_path, metrics_path, *summary_jsons] if p.exists()]
        return {"path_list": files, "has_all": score_path.exists() and metrics_path.exists()}

    def get_deletable(self, chunk, name):
        return []


def _safe_pair_value(matrix):
    """Return matrix[0][1] or None if the slot is missing/non-numeric."""
    if not matrix or len(matrix) < 1 or len(matrix[0]) < 2:
        return None
    v = matrix[0][1]
    return float(v) if isinstance(v, (int, float)) else None
