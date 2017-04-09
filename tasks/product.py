from functools import reduce
from utils import submission as sub
from utils import meta
from tqdm import tqdm  # type: ignore
import pandas as pd  # type: ignore
import numpy as np  # type: ignore

TARGETS = list(sub.MAP_TARGET.keys())
REGIONS = list(sub.MAP_REGION.keys())


def generate_output(kcde: sub.Submission,
                    kde: sub.Submission,
                    sarima: sub.Submission) -> sub.Submission:
    sections = []
    for target in TARGETS:
        for region in REGIONS:
            kcde_point, kcde_bins = kcde.get_subset(region, target)
            kde_point, kde_bins = kde.get_subset(region, target)
            sarima_point, sarima_bins = sarima.get_subset(region, target)

            # Take product of the distributions
            out_bin_X = reduce(np.multiply, [
                kcde_bins.iloc[:, -1].values, kde_bins.iloc[:, -1].values,
                sarima_bins.iloc[:, -1].values
            ])

            # Normalizing
            out_bin_X = out_bin_X / out_bin_X.sum()

            if target in [
                    "one_week", "two_weeks", "three_weeks", "four_weeks",
                    "peak"
            ]:
                # Skip the last bin (13-100%) which can be large
                point_index = np.argmax(out_bin_X[:-1])
            else:
                point_index = np.argmax(out_bin_X)

            point = kcde_bins.iloc[point_index, 0]

            sections.append(sub.build_sub_df(out_bin_X, point, region, target))

    return sub.Submission(pd.concat(sections))


for i in tqdm(range(len(snakemake.input.KCDE))):
    kcde_sub = sub.Submission(snakemake.input.KCDE[i])
    kde_sub = sub.Submission(snakemake.input.KDE[i])
    sarima_sub = sub.Submission(snakemake.input.SARIMA[i])

    out = generate_output(kcde_sub, kde_sub, sarima_sub)
    out.to_csv(snakemake.output.files[i])

# Write metadata file
metadata = {
    "name": "Product of experts model",
    "description": "Product of distributions from KCDE, KDE and SARIMA",
    "url": "#"
}

meta.write_meta(metadata, snakemake.output.meta)
