from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from zipfile import ZipFile
import subprocess
import shutil
import os

HTTP = HTTPRemoteProvider()

SEASON_ROOT = "data/2016-2017"

# Assuming KCDE files are complete
weeks = [p.stem for p in Path(SEASON_ROOT + "/KCDE").glob("*.csv")]

# Inputs
KCDE_files = SEASON_ROOT + "/KCDE/{week}.csv"
KDE_files = SEASON_ROOT + "/KDE/{week}.csv"
SARIMA_files = SEASON_ROOT + "/SARIMA/{week}.csv"

# Outputs
Average_files = SEASON_ROOT + "/Average/{week}.csv"
Average_meta = SEASON_ROOT + "/Average/meta.yml"
NN_files = SEASON_ROOT + "/NN/{week}.csv"
NN_meta = SEASON_ROOT + "/NN/meta.yml"

rule nn:
    input:
        KCDE = expand(KCDE_files, week=weeks),
        KDE = expand(KDE_files, week=weeks),
        SARIMA = expand(SARIMA_files, week=weeks)
    output:
        files = expand(NN_files, week=weeks),
        meta = NN_meta
    message: "Running neural network model"
    script: "tasks/nn.py"

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

rule pull_data:
    input:
        HTTP.remote("github.com/reichlab/flusight/archive/master.zip", keep_local=False, allow_redirects=True)
    output:
        KCDE = SEASON_ROOT + "/KCDE",
        KDE = SEASON_ROOT + "/KDE",
        SARIMA = SEASON_ROOT + "/SARIMA"
    message: "Pulling in latest submission files from flusight"
    run:
        try:
            shutil.rmtree(output.KCDE)
            shutil.rmtree(output.KDE)
            shutil.rmtree(output.SARIMA)
        except:
            pass

        with ZipFile(input[0]) as zf:
            zf.extractall()

        downloaded_data_root = "flusight-master/" + SEASON_ROOT
        shutil.copytree(downloaded_data_root + "/KCDE", output.KCDE)
        shutil.copytree(downloaded_data_root + "/KDE", output.KDE)
        shutil.copytree(downloaded_data_root + "/SARIMA", output.SARIMA)

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
