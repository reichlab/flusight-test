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

            # TODO Check for probability summation
            out_bin_X = np.mean(
                [
                    kcde_bins.iloc[:, -1], kde_bins.iloc[:, -1],
                    sarima_bins.iloc[:, -1]
                ],
                axis=0)

            # TODO Create custom point prediction
            sections.append(
                sub.build_sub_df(out_bin_X, kcde_point, region, target))

    return sub.Submission(pd.concat(sections))


for i in tqdm(range(len(snakemake.input.KCDE))):
    kcde_sub = sub.Submission(snakemake.input.KCDE[i])
    kde_sub = sub.Submission(snakemake.input.KDE[i])
    sarima_sub = sub.Submission(snakemake.input.SARIMA[i])

    out = generate_output(kcde_sub, kde_sub, sarima_sub)
    out.to_csv(snakemake.output.files[i])

# Write metadata file
metadata = {
    "name": "Averaging model",
    "description": "Average predictions from KCDE, KDE and SARIMA",
    "url": "#"
}

meta.write_meta(metadata, snakemake.output.meta)
