__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from tempfile import gettempdir, TemporaryDirectory

from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fastas = snakemake.input.get("fastas")
assert fastas is not None, "input: fastas is a required input parameter"

fastas = [fastas] if isinstance(fastas, str) else fastas
fastas = [f"'{f}'" for f in fastas]

gtf = snakemake.input.get("gtf", "")
if gtf:
    assert gtf.endswith((".gtf", ".gff3")), "input: gtf extension not .gtf/gff3"
    gtf = f"--sjdbGTFfile {gtf}"
    if gtf.endswith(".gff3"):
        gtf += " --sjdbGTFtagExonParentTranscript Parent"

extra = snakemake.params.get("extra", "")

makedirs(snakemake.output[0])

tmp_base_dir = snakemake.resources.get("tmp_dir", gettempdir())

with TemporaryDirectory(dir=tmp_base_dir) as tmp_dir:
    shell(
        "STAR"
        " --runThreadN {snakemake.threads}"
        " --runMode genomeGenerate"
        " --genomeFastaFiles {fastas}"
        " --genomeDir {snakemake.output}"
        " --outTmpDir {tmp_dir}/tmp"
        " {gtf}"
        " {extra}"
        " {log}"
    )
