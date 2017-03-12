from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from zipfile import ZipFile
import subprocess
import shutil
import os

HTTP = HTTPRemoteProvider()

# Assuming KCDE files are complete
weeks = [p.stem for p in Path("data/2016-2017/KCDE").glob("*.csv")]

# Inputs
KCDE_files = "data/2016-2017/KCDE/{week}.csv"
KDE_files = "data/2016-2017/KDE/{week}.csv"
SARIMA_files = "data/2016-2017/SARIMA/{week}.csv"

# Outputs
Average_files = "data/2016-2017/Average/{week}.csv"
Average_meta = "data/2016-2017/Average/meta.yml"

rule average:
    input:
        KCDE = expand(KCDE_files, week=weeks),
        KDE = expand(KDE_files, week=weeks),
        SARIMA = expand(SARIMA_files, week=weeks)
    output:
        files = expand(Average_files, week=weeks),
        meta = Average_meta
    message: "Running averaging model"
    script: "tasks/average.py"

rule all:
    input:
        rules.average.output

rule flusight:
    input:
        HTTP.remote("github.com/reichlab/flusight/archive/master.zip", keep_local=False, allow_redirects=True)
    output: "flusight-local"
    message: "Setting up flusight"
    run:
        try:
            shutil.rmtree(output[0])
        except:
            pass

        with ZipFile(input[0]) as zf:
            zf.extractall()
        os.rename("flusight-master", output[0])

        shutil.rmtree(os.path.join(output[0], "data"))
        os.remove(os.path.join(output[0], "config.yaml"))

        os.symlink("../data", os.path.join(output[0], "data"), target_is_directory=True)
        os.symlink("../config.yaml", os.path.join(output[0], "config.yaml"))

        # Install deps
        os.chdir(output[0])
        subprocess.run(["yarn"])
        os.chdir("..")
