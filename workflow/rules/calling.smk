rule prosolo_calling:
    input:
        single_cell="results/recal/{sc}.sorted.bam",
        single_cell_index="results/recal/{sc}.sorted.bai",
        bulk="results/recal/{individual}.merged_bulk.sorted.bam",
        bulk_index="results/recal/{individual}.merged_bulk.sorted.bai",
        ref=rules.create_full_reference.output,
        ref_idx=rules.full_reference_faidx.output,
        candidates="results/candidate-calls/{individual}.freebayes.{scatteritem}.bcf",
    output:
        "results/calls/{individual}/{sc}.merged_bulk.prosolo.{scatteritem}.bcf",
    params:
        extra="",
    threads: 1
    log:
        "logs/prosolo/{individual}/{sc}.merged_bulk.prosolo.{scatteritem}.log",
    wrapper:
        "v1.3.1/bio/prosolo/single-cell-bulk"


rule sort_calls:
    input:
        "results/calls/{individual}/{sc}.merged_bulk.prosolo.{scatteritem}.bcf",
    output:
        temp("results/calls/{individual}/{sc}.merged_bulk.prosolo.{scatteritem}.sorted.bcf"),
    log:
        "logs/bcf-sort/{individual}/{sc}.merged_bulk.prosolo.{scatteritem}.sorted.log",
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb=8000,
    shell:
        "bcftools sort --max-mem {resources.mem_mb}M --temp-dir `mktemp -d` "
        "-Ob {input} > {output} 2> {log}"


rule bcftools_index_scattered_calls:
    input:
        "results/calls/{individual}/{sc}.merged_bulk.prosolo.{scatteritem}.sorted.bcf",
    output:
        temp("results/calls/{individual}/{sc}.merged_bulk.prosolo.{scatteritem}.sorted.bcf.csi"),
    log:
        "logs/bcftools_index/{individual}/{sc}.merged_bulk.prosolo.{scatteritem}.sorted.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index {input} 2> {log}"


rule gather_scattered_calls:
    input:
        calls=gather.calling("results/calls/{{individual}}/{{sc}}.merged_bulk.prosolo.{scatteritem}.sorted.bcf"),
        indexes=gather.calling("results/calls/{{individual}}/{{sc}}.merged_bulk.prosolo.{scatteritem}.sorted.bcf.csi"),
    output:
        "results/final-calls/{individual}/{sc}.merged_bulk.prosolo.sorted.bcf",
    log:
        "logs/gather_scattered_calls/{individual}/{sc}.merged_bulk.prosolo.sorted.log",
    params:
        extra="-a",
    wrapper:
        "v1.3.1/bio/bcftools/concat"
