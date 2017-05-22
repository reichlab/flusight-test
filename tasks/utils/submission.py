"""
Module to work with submission files
"""

from typing import Dict, List

import numpy as np # type: ignore
import pandas as pd # type: ignore

SUB_HEADER = [
    "Location",
    "Target",
    "Type",
    "Unit",
    "Bin_start_incl",
    "Bin_end_notincl",
    "Value"
]

# Map from targets used in models to that in submissions
MAP_TARGET = {
    model: submission
    for model, submission in zip([
        "onset", "peak_week", "peak", "one_week", "two_weeks",
        "three_weeks", "four_weeks"
    ], [
        "Season onset", "Season peak week", "Season peak percentage",
        "1 wk ahead", "2 wk ahead", "3 wk ahead", "4 wk ahead"
    ])
}

# Map from region used in models to that in submissions
MAP_REGION = {
    "".join(submission.split()[1:]): submission
    for submission in [
            "US National", "HHS Region 1", "HHS Region 2", "HHS Region 3",
            "HHS Region 4", "HHS Region 5", "HHS Region 6", "HHS Region 7",
            "HHS Region 8", "HHS Region 9", "HHS Region 10"
    ]
}


def build_sub_df(bins: np.ndarray, point_prediction, region: str, target: str) -> Dict:
    """
    Create a segment of rows going into submission
    """

    df = {key: [] for key in SUB_HEADER}

    def _append_row(row):
        df[SUB_HEADER[0]].append(MAP_REGION[row[0]])
        df[SUB_HEADER[1]].append(MAP_TARGET[row[1]])
        for i in range(2, 7):
            df[SUB_HEADER[i]].append(row[i])

    if target in ["onset", "peak_week"]:
        # It should have 33 rows
        # TODO Fix the number of bins if the season has different number of weeks
        if bins.shape != (33,):
            raise ValueError(f"bins.shape needed (33,), got {bins.shape}")

        # Not outputting point predictions since that is autocalculated by
        # flusight
        # TODO
        _append_row([region, target, "Point", "week", None, None, point_prediction])

        bin_starts = list(range(40, 53)) + list(range(1, 21))
        bin_ends: List[float] = [i + 1 for i in bin_starts]

        for idx, x in enumerate(bins):
            _append_row([
                region, target, "Bin", "week",
                str(bin_starts[idx]), str(bin_ends[idx]), x
            ])

        # Add none bin with 0 probability
        # TODO
        if target == "onset":
            _append_row([region, target, "Bin", "week", "none", "none", 0.0])

    else:
        # Assume percentage bin of size 131
        if bins.shape != (131,):
            raise ValueError(f"bins.shape needed (131,), got {bins.shape}")

        # Point prediction
        # TODO
        _append_row([region, target, "Point", "percent", None, None, point_prediction])

        bin_starts = np.linspace(0, 13, 131)
        bin_ends = [i + 0.1 for i in bin_starts]
        # Set last bin_end to 100 to conform with submission format
        bin_ends[-1] = 100

        for idx, x in enumerate(bins):
            _append_row([
                region, target, "Bin", "percent",
                f"{bin_starts[idx]:.1f}", f"{bin_ends[idx]:.1f}", x
            ])

    return pd.DataFrame(df)[SUB_HEADER]


class Submission:
    """
    Class for submission file in long format
    """

    def __init__(self, source) -> None:
        """
        Create submission object from df
        """

        try:
            self.df = pd.read_csv(source)
        except Exception:
            self.df = source

    def to_csv(self, file_name) -> None:
        """
        Write submission file to given path
        """

        self.df.to_csv(file_name, na_rep="NA", index=False)

    def get_subset(self, region: str, target: str):
        """
        Return bins and point value for asked region and target
        """

        subset = self.df[
            (self.df["Location"] == MAP_REGION[region]) & \
            (self.df["Target"] == MAP_TARGET[target])
        ]

        bins = subset[["Bin_start_incl", "Bin_end_notincl", "Value"]].iloc[1:, :]
        point = subset["Value"].iloc[0]

        if target == "onset":
            # Don't return none bin
            bins = bins.iloc[:-1, :]

        return point, bins
