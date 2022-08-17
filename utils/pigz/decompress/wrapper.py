__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

file_type = (
    "--zip"
    if snakemake.input[0].endswith((".zip", ".ZIP"))
    else "--zlib"
    if snakemake.input[0].endswith((".zz", ".ZZ"))
    else ""
)

shell(
    "pigz"
    " --processes {snakemake.threads}"
    " {file_type}"
    " --decompress"
    " --stdout {snakemake.input[0]} > {snakemake.output[0]}"
    " {log}"
)
