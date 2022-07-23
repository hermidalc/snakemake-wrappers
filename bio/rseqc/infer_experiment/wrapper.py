__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

import numpy as np
import re
from snakemake.shell import shell

stranded_p_thres = 90
unstranded_p_thres = 55

log = snakemake.log_fmt_shell(stderr=True, stdout=False)

bed = snakemake.input.get("bed")
assert bed is not None, "input: bed is a required input parameter"
bam = snakemake.input.get("bam")
assert bam is not None, "input: bam is a required input parameter"

infer_file = snakemake.output.get("infer")
assert infer_file is not None, "output: infer is a required output parameter"
strand_file = snakemake.output.get("strand")
assert strand_file is not None, "output: strand is a required output parameter"

sample_size = snakemake.params.get("sample_size", 200000)

shell("infer_experiment.py -r {bed} -i {bam} -s {sample_size} > {infer_file} {log}")

# Example:
#
# This is PairEnd Data
# Fraction of reads failed to determine: 0.0234
# Fraction of reads explained by "1++,1--,2+-,2-+": 0.0028
# Fraction of reads explained by "1+-,1-+,2++,2--": 0.9738

with open(snakemake.log[0], "w") as f_log:
    f_log.write(f"Getting strandedness from {input[0]}\n")
    fwd_str = "1++,1--,2+-,2-+"
    rev_str = "1+-,1-+,2++,2--"
    fwd_str_esc = re.escape(fwd_str)
    rev_str_esc = re.escape(rev_str)
    frac_line_regex = re.compile(
        f'^Fraction of reads explained by "({fwd_str_esc}|{rev_str_esc})": (0\.[0-9]+)$',
        flags=re.IGNORECASE,
    )
    fracs = []
    with open("test.txt", "r") as f_in:
        for line in f_in:
            if m := re.match(frac_line_regex, line.strip()):
                fracs.append(float(m.group(2)))
    num_fracs = len(fracs)
    assert num_fracs == 2, f"Output file has {num_fracs} but should have 2"

    percents = np.rint(np.array(fracs) * 100)

    if np.all(percents < unstranded_p_thres):
        stranded = 0
    elif np.any(percents > stranded_p_thres) and np.any(
        percents < 100 - stranded_p_thres
    ):
        stranded = np.where(percents > stranded_p_thres)[0][0] + 1
    else:
        raise ValueError("Cannot determine strandedness from output file")

    with open(strand_file, "w") as f_out:
        f_out.write(
            "forward" if stranded == 1 else "reverse" if stranded == 2 else "no"
        )
