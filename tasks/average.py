from utils import submission as sub
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
            kcde_X = kcde.get_X(region, target)
            kde_X = kde.get_X(region, target)
            sarima_X = sarima.get_X(region, target)

            out_X = np.mean([kcde_X, kde_X, sarima_X], axis=0)
            sections.append(sub.build_sub_df(out_X, 0, region, target))

    return sub.Submission(pd.concat(sections))


for i in tqdm(range(len(snakemake.input.KCDE))):
    kcde_sub = sub.Submission(snakemake.input.KCDE[i])
    kde_sub = sub.Submission(snakemake.input.KDE[i])
    sarima_sub = sub.Submission(snakemake.input.SARIMA[i])

    out = generate_output(kcde_sub, kde_sub, sarima_sub)
    out.to_csv(snakemake.output[i])
