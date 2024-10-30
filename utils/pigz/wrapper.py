__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

decompress = snakemake.params.get("decompress", False)
if decompress:
    decompress = "--decompress"

file_type = (
    "--zip"
    if snakemake.input[0].endswith((".zip", ".ZIP"))
    else "--zlib" if snakemake.input[0].endswith((".zz", ".ZZ")) else ""
)

extra = snakemake.params.get("extra", "")

shell(
    "pigz"
    " --processes {snakemake.threads}"
    " {decompress}"
    " {file_type}"
    " {extra}"
    " --stdout {snakemake.input[0]} > {snakemake.output[0]}"
    " {log}"
)
