import csv
import os
import re
from typing import Dict, Tuple, Union, List
from glob import glob
import logging
import operator

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyteomics import mgf
import seaborn as sns
import spectrum_utils.spectrum as sus


class PeptideSpectrumMatch:
    def __init__(
        self,
        scan: Union[int, None] = None,
        sequence: Union[str, None] = None,
        modifications: Dict = dict(),
        retention_time: Union[float, None] = None,
        q_value: Union[float, None] = None,
        score: Union[float, None] = None,
    ):
        self.scan = scan
        self.sequence = sequence
        self.modifications = modifications
        self.retention_time = retention_time
        self.q_value = q_value
        self.score = score

    def __repr__(self):
        return f"PeptideSpectrumMatch(scan={self.scan}, \
sequence={self.sequence}, modifications={self.modifications}, \
retention_time={self.retention_time}, q_value={self.q_value}, \
score={self.score})"

    def add_pout_modified_sequence(
        self, modified_sequence: str, mod_mapping: Union[Dict, None] = None
    ):
        # Remove leading and trailing amino acids
        if "." in modified_sequence:
            modified_sequence = modified_sequence.split(".")[1]

        # Get all mod labels, locations in original string and label lengths
        mod_pattern = r"\[([^]]*)\]"
        mod_strings = re.findall(mod_pattern, modified_sequence)
        if mod_strings:
            mod_label_locations, mod_label_lengths = zip(
                *[
                    (m.start(0), m.end(0) - m.start(0))
                    for m in re.finditer(mod_pattern, modified_sequence)
                ]
            )

            # Substract length of modification labels to get actual mod locations
            mod_locations = np.array(mod_label_locations) - np.cumsum(
                [0] + list(mod_label_lengths)[:-1]
            )

            if mod_mapping:
                mod_strings = [mod_mapping[ms] for ms in mod_strings]
        else:
            mod_locations, mod_strings = [], []

        # Construct PEPREC modifications string
        self.modifications = "|".join(
            ["|".join([str(ml), ms]) for ml, ms in zip(mod_locations, mod_strings)]
        )
        self.sequence = re.sub(mod_pattern, "", modified_sequence)


class Run:
    def __init__(
        self, run_name: Union[str, None] = None, mgf_dir: str = "", pout_dir: str = "",
    ):
        self.run_name = run_name
        self.mgf_dir = mgf_dir
        self.pout_dir = pout_dir
        self.peptide_spectrum_matches = dict()

    def get_pout_filename(self) -> str:
        """Return pout filename based on pout_dir and run_name."""
        return os.path.join(self.pout_dir, self.run_name + ".pout")

    def get_mgf_filename(self) -> str:
        """Return mgf filename based on mgf_dir and run_name."""
        return os.path.join(self.mgf_dir, self.run_name + ".mgf")

    def num_psms(self) -> int:
        """Get number of PSMs in run."""
        return len(self.peptide_spectrum_matches)

    def read_mgf(
        self, mgf_filename: Union[str, None] = None, no_new_psms: bool = False
    ):
        """Read retention times from MGF file."""
        if not mgf_filename:
            mgf_filename = self.get_mgf_filename()
        with mgf.MGF(mgf_filename) as reader:
            for spectrum in reader:
                scan = int(spectrum["params"]["scans"])
                retention_time = float(spectrum["params"]["rtinseconds"])
                if scan not in self.peptide_spectrum_matches:
                    if not no_new_psms:
                        psm = PeptideSpectrumMatch(
                            scan=scan, retention_time=retention_time
                        )
                        self.peptide_spectrum_matches[scan] = psm
                else:
                    self.peptide_spectrum_matches[scan].retention_time = retention_time

    def read_pout(
        self,
        pout_filename: Union[str, None] = None,
        search_engine: str = "msgfplus",
        no_new_psms: bool = False,
        mod_mapping: Union[Dict, None] = None,
    ):
        """Read PSMs from pout file."""
        if not pout_filename:
            pout_filename = self.get_pout_filename()

        with open(pout_filename, "rt") as csv_file:
            reader = csv.DictReader(csv_file, delimiter="\t")
            for line in reader:
                psmid = line["PSMId"]
                if search_engine == "msgfplus":
                    scan = int(psmid.split("_")[-3])
                else:
                    raise ValueError("Unsupported search engine (supported: msgfplus)")

                score = line["score"]
                q_value = line["q-value"]
                modified_sequence = line["peptide"]

                if scan not in self.peptide_spectrum_matches:
                    if not no_new_psms:
                        psm = PeptideSpectrumMatch(
                            scan=scan, score=score, q_value=q_value
                        )
                        psm.add_pout_modified_sequence(modified_sequence, mod_mapping)
                        self.peptide_spectrum_matches[scan] = psm
                else:
                    self.peptide_spectrum_matches[scan].score = score
                    self.peptide_spectrum_matches[scan].q_value = q_value
                    self.peptide_spectrum_matches[scan].add_pout_modified_sequence(
                        modified_sequence, mod_mapping
                    )

    def get_fraction_missing_rt(self) -> float:
        missing_rt = 0
        total_psms = 0
        for _, psm in self.peptide_spectrum_matches.items():
            if not psm.retention_time:
                missing_rt += 1
            total_psms += 1

        return missing_rt / total_psms

    def _min_val_key(self, dict_) -> str:
        """Return key with min value in dict."""
        return min(dict_.items(), key=operator.itemgetter(1))[0]

    def get_best_psms(self) -> Dict[Tuple[str, str], int]:
        """Get best PSM for each unique peptide-modification combination."""
        psm_dict = dict()
        for _, psm in self.peptide_spectrum_matches.items():
            if (psm.sequence, psm.modifications) not in psm_dict:
                psm_dict[(psm.sequence, psm.modifications)] = {psm.scan: psm.q_value}
            else:
                psm_dict[(psm.sequence, psm.modifications)][psm.scan] = psm.q_value

        best_psms = {
            self.peptide_spectrum_matches[self._min_val_key(scans)]
            for _, scans in psm_dict.items()
        }

        return best_psms

    # TODO: Needs optimization to read PSMs and immediatly go to dataframe, skipping
    # PeptideSpectrumMatch objects
    def to_dataframe(self) -> pd.DataFrame:
        """Dump all PSMs into a pandas.DataFrame."""

        idx = range(self.num_psms())
        scans = pd.Series(name="scan", index=idx, dtype=np.uint32)
        sequences = pd.Series(name="sequence", index=idx, dtype=object)
        modifications = pd.Series(name="modifications", index=idx, dtype=object)
        q_values = pd.Series(name="q_value", index=idx, dtype=np.float32)
        scores = pd.Series(name="score", index=idx, dtype=np.float32)
        retention_times = pd.Series(name="retention_time", index=idx, dtype=np.float32)

        for idx, (_, psm) in enumerate(self.peptide_spectrum_matches.items()):
            scans[idx] = psm.scan
            sequences[idx] = psm.sequence
            modifications[idx] = psm.modifications
            retention_times[idx] = psm.retention_time
            q_values[idx] = psm.q_value
            scores[idx] = psm.score

        df = pd.concat(
            [scans, sequences, modifications, retention_times, q_values, scores], axis=1
        )

        return df


class RunCollection:
    def __init__(
        self,
        name: str,
        root_dir: str = "",
        mgf_subdir: str = "mgf",
        pout_subdir: str = "pout",
    ):
        self.name = name
        self.root_dir = root_dir
        self.mgf_subdir = mgf_subdir
        self.pout_subdir = pout_subdir

        self.runs = dict()

    def _add_runs(
        self,
        run_list: List[str],
        read_psms: bool,
        mod_mapping: Union[Dict, None] = None,
    ):
        """Add runs from list of run names."""
        for run in run_list:
            self.runs[run] = Run(
                run_name=run,
                mgf_dir=os.path.join(self.root_dir, self.mgf_subdir),
                pout_dir=os.path.join(self.root_dir, self.pout_subdir),
            )
            if read_psms:
                logging.debug("Reading PSMs for %s", run)
                self.runs[run].read_pout(mod_mapping=mod_mapping)
                self.runs[run].read_mgf(no_new_psms=True)

    def add_runs_by_glob(
        self,
        name_pattern: str = "*",
        read_psms: bool = True,
        mod_mapping: Union[Dict, None] = None,
    ):
        """Add runs by using glob to find all mgf files."""
        mgf_pattern = os.path.join(
            self.root_dir, self.mgf_subdir, name_pattern + ".mgf"
        )
        run_list = glob(mgf_pattern)
        run_list = [os.path.splitext(os.path.basename(run))[0] for run in run_list]
        self._add_runs(run_list, read_psms, mod_mapping=mod_mapping)

    def add_runs_by_list(
        self,
        run_list: List[str],
        read_psms: bool = True,
        mod_mapping: Union[Dict, None] = None,
    ):
        """Add runs from a list of run names."""
        self._add_runs(run_list, read_psms, mod_mapping=mod_mapping)

    def to_dataframe(self) -> pd.DataFrame:
        """Dump all PSMs into a pandas.DataFrame."""
        dfs = []
        for run_name, run in self.runs.items():
            run_df = run.to_dataframe()
            run_df["run"] = run_name
            dfs.append(run_df)
        df = pd.concat(dfs, axis=0)
        df["collection"] = self.name
        return df

    def _calibrate_retention_times(
        self,
        original: List[float],
        original_shared: List[float],
        reference_shared: List[float],
    ) -> List[float]:
        """
        Calibrate retention times to a reference.
        
        Given a list of shared retention times between the orignal set and a reference
        set, calibrate a full list of retention times.
        """
        original = np.array(original, dtype=np.float64)
        original_shared = np.array(original_shared, dtype=np.float64)
        reference_shared = np.array(reference_shared, dtype=np.float64)

        # Insert zero at beginning of set and medians
        original_shared = np.insert(original_shared, 0, 0)
        reference_shared = np.insert(reference_shared, 0, 0)

        # Hard copy all, to apply calibration to
        original_calibrated = original.copy()

        # Calculate differences
        diff = original_shared - reference_shared

        # Loop over all median differences to apply to all
        for i in range(1, len(reference_shared)):
            original_calibrated[
                np.logical_and(
                    original > original_shared[i - 1], original <= original_shared[i]
                )
            ] += diff[i]

        # Also apply last difference to elements higher than the last shared rt
        original_calibrated[original > original_shared[-1]] += diff[-1]

        return original_calibrated

    def calibrate_collection(
        self,
        psms: pd.DataFrame = None,
        top_n: Union[float, None] = None,
        q_value_threshold: float = 0.01,
        plot: bool = False,
    ) -> pd.DataFrame:
        """Calibrate retention times in a collection to one run in the collection."""
        if psms is None:
            psms = self.to_dataframe()

        psms["modifications"].fillna("", inplace=True)
        psms = psms[psms["q_value"] <= q_value_threshold]

        print("calibrating ", self.name)
        print("#PSMs @0.01FDR: ", len(psms))

        # Collapse PSMs to unique sequence/modifications per run with median rt
        gb_object = psms.groupby(["sequence", "modifications", "run", "collection"])
        rt = gb_object["retention_time"].median().rename("retention_time_median")
        q_value = gb_object["q_value"].mean().rename("q_value_mean")
        psms_medians = pd.concat([rt, q_value], axis=1).reset_index()

        print("#Peptidoforms: ", len(psms_medians))
        # Get number of runs in which a peptide-mod is
        run_counts = (
            psms_medians.groupby(["sequence", "modifications", "collection"])
            .size()
            .rename("run_counts")
            .reset_index()
        )
        psms_medians = psms_medians.merge(
            run_counts, on=["sequence", "modifications", "collection"]
        )

        # Filter on shared PSMs only
        psms_medians_shared = psms_medians[
            psms_medians["run_counts"] == len(psms_medians["run"].unique())
        ]

        if len(psms_medians_shared) > 0:
            ref_run = psms_medians_shared["run"].iloc[0]
            print("Reference run: ", ref_run)

            if top_n:
                psms_medians_shared = psms_medians_shared.sort_values(
                    "q_value_mean", ascending=True
                ).head(top_n)

            psms_medians_shared = psms_medians_shared[
                psms_medians_shared["run"] == ref_run
            ][["sequence", "modifications", "retention_time_median"]]
            psms_medians_shared = psms_medians_shared.rename(
                columns={"retention_time_median": "retention_time_reference"}
            )
            psms_medians = psms_medians.merge(psms_medians_shared, how="left")

            print("#Shared peptidoforms: ", len(psms_medians_shared))

            calibrated = []

            for run in psms_medians["run"].unique():
                psms_run = psms_medians[psms_medians["run"] == run].copy()

                # Get arrays with original, original_shared, and reference_shared
                psms_run_shared = psms_run[
                    ~psms_run["retention_time_reference"].isna()
                ].sort_values("retention_time_median")
                original_shared = psms_run_shared["retention_time_median"]
                reference_shared = psms_run_shared["retention_time_reference"]
                original = psms_run["retention_time_median"].sort_values()

                # Sort by retention_time_median, before adding calibrated rt's!!!
                psms_run = psms_run.sort_values("retention_time_median")

                # calibrate
                psms_run["retention_time_calibrated"] = self._calibrate_retention_times(
                    original, original_shared, reference_shared
                )

                calibrated.append(psms_run)

            psms_calibrated = pd.concat(calibrated, axis=0)

            if plot:
                plt.figure()
                sns.lmplot(
                    data=psms_calibrated,
                    x="retention_time_calibrated",
                    y="retention_time_reference",
                    hue="run",
                    scatter_kws={"s": 4},
                    fit_reg=False,
                )
                plt.savefig(self.name + ".png")

        else:
            psms_calibrated = psms_medians
            psms_calibrated["retention_time_calibrated"] = psms_calibrated[
                "retention_time_median"
            ]

        # Calculate medians of calibrated retention times
        psms_calibrated_medians = (
            psms_calibrated.groupby(["sequence", "modifications"])[
                "retention_time_calibrated"
            ]
            .median()
            .reset_index()
        )

        return psms_calibrated_medians


# TODO: Implement CLI
def main():
  pass
  
if __name__ == "__main__":
    main()
