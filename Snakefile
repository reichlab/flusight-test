from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

# Assuming KCDE files are complete
weeks = [p.stem for p in Path("data/2016-2017/KCDE").glob("*.csv")]

# Inputs
KCDE_files = "data/2016-2017/KCDE/{week}.csv"
KDE_files = "data/2016-2017/KDE/{week}.csv"
SARIMA_files = "data/2016-2017/SARIMA/{week}.csv"

# Outputs
Average_files = "data/2016-2017/Average/{week}.csv"

rule average:
    input:
        KCDE = expand(KCDE_files, week=weeks),
        KDE = expand(KDE_files, week=weeks),
        SARIMA = expand(SARIMA_files, week=weeks)
    output: expand(Average_files, week=weeks)
    message: "Running averaging model"
    script: "tasks/average.py"

rule all:
    input:
        rules.average.output

rule flusight:
    input:
        HTTP.remote("github.com/reichlab/flusight/archive/master.zip", keep_local=False, allow_redirects=True)
    output: touch("flusight.done")
    message: "Setting up flusight"
    run: ...
