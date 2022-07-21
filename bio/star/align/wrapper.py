__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

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
    fq_str = f"'{fq1_str} {fq2_str}'"
else:
    fq_str = f"'{fq1_str}'"


if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""

index = snakemake.input.get("index", snakemake.params.get("index", "GenomeDir/"))

gtf = snakemake.input.get("gtf")
assert gtf is not None, "input: gtf is a required input parameter"
assert gtf.endswith(".gtf"), "input: gtf extension not .gtf"

readlength = snakemake.input.get("readlength")
sjdb_overhang = readlength - 1 if readlength else "100"

extra = snakemake.params.get("extra", "")
sjdb = snakemake.input.get("sjdb")
if sjdb:
    extra += f" --sjdbFileChrStartEnd '{sjdb}'"

with TemporaryDirectory() as tmp_dir:
    shell(
        "STAR"
        " --runThreadN {snakemake.threads}"
        " --readFilesIn {fq_str}"
        " --genomeDir {index}"
        " --outFileNamePrefix {tmp_dir}/"
        " --outTmpDir {tmp_dir}/tmp"
        " --sjdbGTFfile {gtf}"
        " --sjdbOverhang {sjdb_overhang}"
        " {readcmd}"
        " {extra}"
        " {log}"
    )

    if "SortedByCoordinate" in extra:
        bam_prefix = "Aligned.sortedByCoord.out"
    else:
        bam_prefix = "Aligned.out"

    if snakemake.output.get("bam"):
        shell("cat {tmpdir}/{bam_prefix}.bam > {snakemake.output.bam:q}")
    if snakemake.output.get("sam"):
        shell("cat {tmpdir}/{bam_prefix}.sam > {snakemake.output.sam:q}")
    if snakemake.output.get("read_counts"):
        shell("cat {tmpdir}/ReadsPerGene.out.tab > {snakemake.output.reads_per_gene:q}")
    if snakemake.output.get("chim_junc"):
        shell("cat {tmpdir}/Chimeric.out.junction > {snakemake.output.chim_junc:q}")
    if snakemake.output.get("sj"):
        shell("cat {tmpdir}/SJ.out.tab > {snakemake.output.sj:q}")
    if snakemake.output.get("log"):
        shell("cat {tmpdir}/Log.out > {snakemake.output.log:q}")
    if snakemake.output.get("log_prog"):
        shell("cat {tmpdir}/Log.progress.out > {snakemake.output.log_progress:q}")
    if snakemake.output.get("log_final"):
        shell("cat {tmpdir}/Log.final.out > {snakemake.output.log_final:q}")
