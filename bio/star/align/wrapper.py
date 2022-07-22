__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

import re
from os.path import join
from shutil import rmtree
from tempfile import TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fq1 = snakemake.input.get("fq1")
assert fq1 is not None, "input: fq1 is a required input parameter"
fq1 = (
    [snakemake.input.fq1]
    if isinstance(snakemake.input.fq1, str)
    else snakemake.input.fq1
)
fq1_str = ",".join(fq1)

fq2 = snakemake.input.get("fq2")
if fq2:
    fq2 = (
        [snakemake.input.fq2]
        if isinstance(snakemake.input.fq2, str)
        else snakemake.input.fq2
    )
    assert len(fq1) == len(fq2), "input: equal number of files required for fq1 and fq2"
    fq2_str = ",".join(fq2)
    fq_str = f"'{fq1_str}' '{fq2_str}'"
else:
    fq_str = f"'{fq1_str}'"


if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""

index = snakemake.input.get("index", snakemake.params.get("index", "GenomeDir/"))

out_dir = snakemake.params.get("out_dir")
assert out_dir is not None, "params: out_dir is a required parameter"


gtf = snakemake.input.get("gtf")
assert gtf is not None, "input: gtf is a required input parameter"
assert gtf.endswith(".gtf"), "input: gtf extension not .gtf"

readlength = snakemake.params.get("readlength")
if readlength is None:
    readlength_file = snakemake.input.get("readlength")
    assert (
        readlength_file is not None
    ), "input/params: readlength is a required parameter"
    with open(readlength_file, "r") as fh:
        readlength = re.sub("\D+", "", fh.readline())

sjdb_overhang = int(readlength) - 1 if readlength else 100

extra = snakemake.params.get("extra", "")
sj = snakemake.input.get("sj")
if sj:
    extra += f" --sjdbFileChrStartEnd {sj}"

with TemporaryDirectory() as tmp_dir:
    shell(
        "STAR"
        " --runThreadN {snakemake.threads}"
        " --readFilesIn {fq_str}"
        " --genomeDir {index}"
        " --outFileNamePrefix {out_dir}/"
        " --outTmpDir {tmp_dir}/tmp"
        " --sjdbGTFfile {gtf}"
        " --sjdbOverhang {sjdb_overhang}"
        " {readcmd}"
        " {extra}"
        " {log}"
    )

rmtree(join(out_dir, "_STARgenome"), ignore_errors=True)
