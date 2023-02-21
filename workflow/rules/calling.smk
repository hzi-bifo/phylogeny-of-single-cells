rule scatter_candidate_calls:
    input:
        "results/candidate_calls/{individual}.freebayes.norm.bcf",
    output:
        scatter.prosolo_chunks(
            "results/candidate_calls/prosolo_chunks/{{individual}}/{scatteritem}.candidate_calls.bcf"
        ),
    log:
        "logs/candidate_calls/prosolo_chunks/{individual}.scatter.candidate_calls.log",
    conda:
        "../envs/rbt.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 40 - 1,
        mem_mb=lambda wildcards, input, attempt: input.size_mb * 1.4 * attempt,
    shell:
        "rbt vcf-split {input} {output}"


rule prosolo_calling:
    input:
        single_cell="results/recal/{sc}.sorted.bam",
        single_cell_index="results/recal/{sc}.sorted.bai",
        bulk="results/recal/{individual}.merged_bulk.sorted.bam",
        bulk_index="results/recal/{individual}.merged_bulk.sorted.bai",
        ref=rules.create_full_reference.output,
        ref_idx=rules.full_reference_faidx.output,
        candidates="results/candidate_calls/prosolo_chunks/{individual}/{scatteritem}.candidate_calls.bcf",
    output:
        "results/calls/{individual}/{sc}/{scatteritem}.merged_bulk.prosolo.bcf",
    params:
        extra="",
    threads: 1
    resources:
        runtime=lambda wildcards, attempt: 6 * attempt * 60 - 1,
        mem_mb=lambda wildcards, attempt: attempt * 4000,
    log:
        "logs/prosolo/{individual}/{sc}.{scatteritem}.merged_bulk.prosolo.log",
    wrapper:
        "v1.21.1/bio/prosolo/single-cell-bulk"


rule sort_calls:
    input:
        "results/calls/{individual}/{sc}/{scatteritem}.merged_bulk.prosolo.bcf",
    output:
        "results/calls/{individual}/{sc}/{scatteritem}.merged_bulk.prosolo.sorted.bcf"
    log:
        "logs/bcf-sort/{individual}/{sc}/{scatteritem}.merged_bulk.prosolo.sorted.log",
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb=4000,
        runtime=lambda wildcards, attempt: 40 * attempt - 1,
    shell:
        "bcftools sort --max-mem {resources.mem_mb}M --temp-dir `mktemp -d` "
        "-Ob {input} > {output} 2> {log}"


rule bcftools_index_region_calls:
    input:
        "results/calls/{individual}/{sc}/{scatteritem}.merged_bulk.prosolo.sorted.bcf",
    output:
        "results/calls/{individual}/{sc}/{scatteritem}.merged_bulk.prosolo.sorted.bcf.csi"
    log:
        "logs/bcftools_index/{individual}/{sc}/{scatteritem}.merged_bulk.prosolo.sorted.log",
    conda:
        "../envs/bcftools.yaml"
    resources:
        runtime=lambda wildcards, attempt: 40 * attempt - 1,
    shell:
        "bcftools index {input} 2> {log}"


rule aggregate_prosolo_chunk_calls:
    input:
        calls=gather.prosolo_chunks(
            "results/calls/{{individual}}/{{sc}}/{scatteritem}.merged_bulk.prosolo.sorted.bcf"
        ),
        indexes=gather.prosolo_chunks(
            "results/calls/{{individual}}/{{sc}}/{scatteritem}.merged_bulk.prosolo.sorted.bcf.csi"
        ),
    output:
        "results/calls/{individual}/{sc}.merged_bulk.prosolo.sorted.bcf",
    log:
        "logs/gather_scattered_calls/{individual}/{sc}.merged_bulk.prosolo.sorted.log",
    params:
        extra="-a",
    resources:
        runtime=lambda wildcards, attempt: 40 * attempt - 1,
        mem_mb=lambda wildcards, input, attempt: input.size_mb * (1 + 0.5 * attempt)
    wrapper:
        "v1.21.1/bio/bcftools/concat"


# Until we have integration of the prosolo with the current varlociraptor
# version and all of its useful priors, we have to do this rough FDR control
# to avoid useless / misleading genotype likelihoods at sites where a sample
# has no or almost no data. The heterozygous genotype probability
# (PROB_ADO_TO_REF + PROB_ADO_TO_ALT + PROB_HET) often goes up above 0.6.
rule prosolo_control_fdr:
    input:
        "results/calls/{individual}/{sc}.merged_bulk.prosolo.sorted.bcf",
    output:
        "results/genotype_fdr/{individual}/{sc}.merged_bulk.prosolo.sorted.{genotype}.fdr_controlled.bcf",
    params:
        # comma-separated set of events for whose (joint)
        # false discovery rate you want to control
        events=lambda wc: "ADO_TO_ALT,ADO_TO_REF,HET"
        if wc.genotype == "het"
        else "HOM_ALT,ERR_REF"
        if wc.genotype == "hom_alt"
        else "HOM_REF,ERR_ALT"
        if wc.genotype == "hom_ref"
        else "ONLY-USE-het-hom_alt-or-hom_ref-for-genotype",
        # false discovery rate to control for
        fdr=0.25,
    log:
        "logs/final_calls/{individual}/{sc}.merged_bulk.prosolo.sorted.{genotype}.fdr_controlled.log",
    resources:
        runtime=lambda wildcards, attempt: 90 * attempt - 1,
        mem_mb=lambda wildcards, attempt: 3000 * attempt
    wrapper:
        "v1.21.1/bio/prosolo/control-fdr"
