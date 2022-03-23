rule prosolo_calling:
    input:
        single_cell = "results/recal/{sc}.sorted.bam",
        single_cell_index = "results/recal/{sc}.sorted.bai",
        bulk = "results/recal/{individual}.merged_bulk.sorted.bam",
        bulk_index = "results/recal/{individual}.merged_bulk.sorted.bai",
        ref = rules.create_full_reference.output,
        ref_idx = rules.full_reference_faidx.output,
        candidates = "results/candidate-calls/{individual}.freebayes.bcf",
    output:
        "results/calls/{individual}/{sc}.merged_bulk.prosolo.bcf"
    params:
        extra = ""
    threads:
        1
    log:
        "logs/prosolo/{individual}/{sc}.prosolo.log"
    wrapper:
        "v1.3.1/bio/prosolo/single-cell-bulk"
