__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

import gffutils

print("Converting GTF to BED", flush=True)

db_file = snakemake.output.get("db")
assert db_file is not None, "output: db is a required output parameter"
bed_file = snakemake.output.get("bed")
assert bed_file is not None, "ouitput: bed is a required output parameter"

db = gffutils.create_db(
    snakemake.input[0],
    dbfn=db_file,
    force=True,
    keep_order=True,
    merge_strategy="merge",
    sort_attribute_values=True,
    disable_infer_genes=True,
    disable_infer_transcripts=True,
)

with open(bed_file, "w") as fh:
    for tx in db.features_of_type("transcript", order_by="start"):
        bed = [s.strip() for s in db.bed12(tx).split("\t")]
        bed[3] = tx.id
        fh.write("{}\n".format("\t".join(bed)))
