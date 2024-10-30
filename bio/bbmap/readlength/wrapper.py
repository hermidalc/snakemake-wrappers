__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

fq1 = snakemake.input.get("fq1") or snakemake.input[0]
assert fq1 is not None, "input: fq1 is a required named or positional input parameter"

fq2 = snakemake.input.get("fq2") or snakemake.input[1]
if fq2:
    in2 = f"in2='{fq2}'"
else:
    in2 = ""

extra = snakemake.params.get("extra", "")

shellcmd = (
    f"JAVA_OPTS='-XX:ActiveProcessorCount={snakemake.threads}' readlength.sh"
    f" in='{fq1}'"
    f" {in2}"
    f" bin=1"
    f" out='{snakemake.output}'"
    f" {extra}"
    f" {log}"
)
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")
shell(shellcmd)
