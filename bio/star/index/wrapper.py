__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from tempfile import TemporaryDirectory
from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")

makedirs(snakemake.output[0])

with TemporaryDirectory() as tmp_dir:
    shell(
        "STAR"
        " --runThreadN {snakemake.threads}"
        " --runMode genomeGenerate"
        " --genomeFastaFiles {snakemake.input}"
        " --genomeDir {snakemake.output}"
        " --outTmpDir {tmp_dir}/tmp"
        " {extra}"
        " {log}"
    )
