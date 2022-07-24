__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

import pandas as pd

count_files = snakemake.input.get("counts")
assert count_files is not None, "input: count_files is a required parameter"
sample_names = snakemake.params.get("samples")
assert sample_names is not None, "params: samples is a required parameter"

strand = snakemake.params.get("strand")
if strand is not None:
    strands = [strand] if isinstance(strand, str) else strand
else:
    strand_files = snakemake.input.get("strand")
    assert strand_files is not None, "input/params: strand is a required parameter"
    strands = []
    for strand_file in strand_files:
        with open(strand_file, "r") as fh:
            strands.append(fh.readline().strip())

count_matrix = pd.DataFrame()
for count_file, sample_name, strand in zip(count_files, sample_names, strands):
    strand_idx = 2 if strand in ("forward", "yes") else 3 if strand == "reverse" else 1
    counts = pd.read_csv(
        count_file, sep="\t", header=None, index_col=0, usecols=[0, strand_idx]
    )
    counts.columns = [sample_name]
    count_matrix = pd.concat([count_matrix, counts], axis=1, verify_integrity=True)
    assert count_matrix.shape[0] == counts.shape[0], "Count files do not have same rows"

count_matrix = count_matrix.loc[~count_matrix.index.str.startswith("N_")]
count_matrix.index.name = "ID_REF"
count_matrix.to_csv(snakemake.output[0], sep="\t")
