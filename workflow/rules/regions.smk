rule build_sample_regions:
    input:
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bai",
    output:
        "results/regions/{individual}/{sample}.mosdepth.global.dist.txt",
        "results/regions/{individual}/{sample}.quantized.bed.gz",
        summary="results/regions/{individual}/{sample}.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth/regions/{individual}_{sample}.log",
    params:
        extra="--no-per-base",
        quantize="1:",
    resources:
        runtime=lambda wc, attempt: attempt * 30 - 1,
    wrapper:
        "v1.21.1/bio/mosdepth"


# TODO: use this rule instead of build_sample_regions
# TODO: make minimum coverage per-sample configurable via config
rule build_sample_regions_cov_2:
    input:
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bai",
    output:
        "results/regions/{individual}/{sample}.cov_2.mosdepth.global.dist.txt",
        "results/regions/{individual}/{sample}.cov_2.quantized.bed.gz",
        summary="results/regions/{individual}/{sample}.cov_2.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth/regions/{individual}_{sample}.cov_2.log",
    params:
        extra="--no-per-base",
        quantize="2:",
    resources:
        runtime=lambda wc, attempt: attempt * 30 - 1,
    wrapper:
        "v1.21.1/bio/mosdepth"



rule merge_individual_regions:
    input:
        lambda wc: expand(
            "results/regions/{{individual}}/{sample}.quantized.bed.gz",
            sample=get_individual_samples(wc.individual),
        ),
    output:
        "results/regions/{individual}.covered_regions.bed",
    log:
        "logs/regions/{individual}_covered_regions.log",
    conda:
        "../envs/bedtools.yaml"
    resources:
        runtime=59,
    shell:
        "zcat {input} | sort -k1,1 -k2,2n - | mergeBed -i - -d 15000 > {output} 2> {log}"


rule filter_individual_regions:
    input:
        covered="results/regions/{individual}.covered_regions.bed",
        fai=rules.full_reference_faidx.output,
    output:
        "results/regions/{individual}.covered_regions.filtered.bed",
    params:
        chroms=config["ref"].get("n_chromosomes", 25),
    log:
        "logs/regions/{individual}.covered_regions.filtered.log",
    shell:
        "cat {input.covered} | grep -f <(head -n {params.chroms} {input.fai} | "
        'awk \'{{print "^"$1"\\t"}}\') '
        "> {output} 2> {log}"
