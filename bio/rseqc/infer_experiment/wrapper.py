__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

bed = snakemake.input.get("bed")
assert bed is not None, "input: bed is a required input parameter"
bam = snakemake.input.get("bam")
assert bam is not None, "input: bam is a required input parameter"

shell("infer_experiment.py -r {bed} -i {bam} > {snakemake.output} {log}")
