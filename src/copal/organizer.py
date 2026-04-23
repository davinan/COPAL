"""Move freshly produced AF3/Boltz2/RF3 outputs into a hierarchical
predictions store, keeping a JSON log of bytes saved/deleted per run.
"""
from __future__ import annotations

import json
import os
import time
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

from copal.prediction_path_managers import AF3PathManager, Boltz2PathManager, RF3PathManager


class Organizer:
    def __init__(self, manager_dict, dst_dir, log_path: Path | None = None):
        self.manager_dict = manager_dict
        self.dst_dir = Path(dst_dir)
        self.log_path = log_path

    def load_log(self):
        if self.log_path is None or not self.log_path.exists():
            return {}
        with open(self.log_path) as fh:
            return json.load(fh)

    def save_log(self, data):
        if self.log_path is None:
            return
        with open(self.log_path, "w") as fh:
            json.dump(data, fh)

    @staticmethod
    def _split_name(name: str) -> tuple[str, str]:
        """Split `<receptor>_<ligand>` where receptor may contain underscores.

        Prefers the last `_AHL` (paper benchmark uses `AHL01`..`AHL58`);
        otherwise falls back to the last underscore so arbitrary
        receptor/ligand IDs still work.
        """
        if "_AHL" in name:
            luxr, ahl = name.rsplit("_AHL", 1)
            return luxr, "AHL" + ahl
        luxr, ahl = name.rsplit("_", 1)
        return luxr, ahl

    def get_range_string(self, name: str, batch_size: int = 200) -> str:
        """Legacy chunk-range prefix used only for `p###`-style receptor IDs.

        Returns "all" for any other receptor naming scheme, producing a
        flat layout `<dst_dir>/<LUXR>/{structures,confidences}/...` that
        the analyzer's `basepath.glob("*")` iterates directly.
        """
        suffix = name[5:]
        if len(name) > 5 and suffix.isdigit():
            batch_id = int(suffix)
            i = batch_id // batch_size
            return f"{batch_size * i}_to_{batch_size * i + batch_size - 1}"
        return "all"

    def organize(self, stage_only: bool = False):
        curr_log = self.load_log()
        time_string = datetime.now().strftime("%H:%M %Y-%m-%d")
        curr_log[time_string] = {}

        for manager_name, manager in self.manager_dict.items():
            curr_log[time_string][manager_name] = {}
            count = 0
            save_sum = 0.0
            delete_sum = 0.0
            t0 = time.time()
            print(f"Organizing: {manager_name}")

            for chunk, name in manager.pred_dir_iter():
                count += 1
                ok, org = manager.get(chunk, name)
                if stage_only:
                    print(f"name: {name} ok: {ok}")
                    continue

                if self.check_moved(manager_name, name):
                    self.delete_files(org["to_delete"])
                    delete_sum += org["to_delete_total_MB"] / 1024

                if ok:
                    save_sum += org["saving_total_MB"] / 1024
                    delete_sum += org["to_delete_total_MB"] / 1024
                    self.move_structures(manager_name, name, org["structures_to_save"])
                    self.move_confidences(manager_name, name, org["confidence_to_save"])
                    if self.check_moved(manager_name, name):
                        self.delete_files(org["to_delete"])

                if count % 100 == 1:
                    curr_log[time_string][manager_name].update({
                        "saved_GB": save_sum,
                        "deleted_GB": delete_sum,
                        "processed": count,
                        "time_elapsed_min": f"{(time.time() - t0) / 60:.2f} minutes",
                    })
                    self.save_log(curr_log)

            curr_log[time_string][manager_name].update({
                "saved_GB": save_sum,
                "deleted_GB": delete_sum,
                "processed": count,
                "time_elapsed_min": f"{(time.time() - t0) / 60:.2f} minutes",
            })
            self.save_log(curr_log)
            print("Organize_status:")
            print(curr_log[time_string][manager_name])
            print("*" * 100)

    def delete_files(self, file_list):
        for file in file_list:
            os.remove(file)

    def check_moved(self, model_name, name):
        luxr_name, _ahl = self._split_name(name)
        chunk_name = self.get_range_string(luxr_name)
        struct_dir = self.dst_dir / chunk_name / luxr_name.upper() / "structures" / name / model_name
        conf_dir = self.dst_dir / chunk_name / luxr_name.upper() / "confidences" / name / model_name
        if not conf_dir.exists() or not struct_dir.exists():
            return False
        return any(struct_dir.iterdir()) and any(conf_dir.iterdir())

    def _move_into(self, model_name, name, file_list, kind):
        luxr_name, _ahl = self._split_name(name)
        chunk_name = self.get_range_string(luxr_name)
        dest = self.dst_dir / chunk_name / luxr_name.upper() / kind / name / model_name
        dest.mkdir(parents=True, exist_ok=True)
        for file in file_list:
            file.rename(dest / file.name)

    def move_structures(self, model_name, name, file_list):
        self._move_into(model_name, name, file_list, "structures")

    def move_confidences(self, model_name, name, file_list):
        self._move_into(model_name, name, file_list, "confidences")


def main():
    parser = ArgumentParser(description="Organize raw model predictions into a hierarchical store.")
    parser.add_argument("--pred-dir", type=Path, required=True,
                        help="Root containing AF3/, Boltz2/, RF3/ subdirs from inference.")
    parser.add_argument("--dst-dir", type=Path, required=True,
                        help="Destination root for organized predictions.")
    parser.add_argument("--log-path", type=Path, default=None,
                        help="Optional JSON log path for run statistics.")
    parser.add_argument("--models", default="boltz2,af3,rf3",
                        help="Comma-separated subset of {boltz2,af3,rf3} to organize.")
    parser.add_argument("--stage-only", action="store_true")
    args = parser.parse_args()

    factory = {
        "af3": lambda p: AF3PathManager(p / "AF3"),
        "boltz2": lambda p: Boltz2PathManager(p / "Boltz2"),
        "rf3": lambda p: RF3PathManager(p / "RF3"),
    }
    selected = [m.strip() for m in args.models.split(",") if m.strip()]
    man_dict = {m: factory[m](args.pred_dir) for m in selected}
    Organizer(man_dict, args.dst_dir, args.log_path).organize(stage_only=args.stage_only)


if __name__ == "__main__":
    main()
