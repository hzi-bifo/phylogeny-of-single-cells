rule prosolo_calling:
    input:
        single_cell="results/recal/{sc}.sorted.bam",
        single_cell_index="results/recal/{sc}.sorted.bai",
        bulk="results/recal/{individual}.merged_bulk.sorted.bam",
        bulk_index="results/recal/{individual}.merged_bulk.sorted.bai",
        ref=rules.create_full_reference.output,
        ref_idx=rules.full_reference_faidx.output,
        candidates="results/candidate-calls/{individual}.{region}.freebayes.norm.bcf",
    output:
        "results/calls/{individual}/{sc}/{region}.merged_bulk.prosolo.bcf",
    params:
        extra="",
    threads: 1
    resources:
        runtime=lambda wildcards, attempt: f"{attempt}:59:00"    
    log:
        "logs/prosolo/{individual}/{sc}.{region}.merged_bulk.prosolo.log",
    wrapper:
        "v1.21.1/bio/prosolo/single-cell-bulk"


rule sort_calls:
    input:
        "results/calls/{individual}/{sc}/{region}.merged_bulk.prosolo.bcf",
    output:
        temp("results/calls/{individual}/{sc}/{region}.merged_bulk.prosolo.sorted.bcf"),
    log:
        "logs/bcf-sort/{individual}/{sc}/{region}.merged_bulk.prosolo.sorted.log",
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb=4000,
    shell:
        "bcftools sort --max-mem {resources.mem_mb}M --temp-dir `mktemp -d` "
        "-Ob {input} > {output} 2> {log}"


rule bcftools_index_region_calls:
    input:
        "results/calls/{individual}/{sc}/{region}.merged_bulk.prosolo.sorted.bcf",
    output:
        temp("results/calls/{individual}/{sc}/{region}.merged_bulk.prosolo.sorted.bcf.csi"),
    log:
        "logs/bcftools_index/{individual}/{sc}/{region}.merged_bulk.prosolo.sorted.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index {input} 2> {log}"


rule aggregate_prosolo_region_calls:
    input:
        calls=aggregate_prosolo_region_calls_input(),
        indexes=aggregate_prosolo_region_calls_input(ext=".bcf.csi"),
    output:
        "results/final-calls/{individual}/{sc}.merged_bulk.prosolo.sorted.bcf",
    log:
        "logs/gather_scattered_calls/{individual}/{sc}.merged_bulk.prosolo.sorted.log",
    params:
        extra="-a",
    wrapper:
        "v1.21.1/bio/bcftools/concat"
