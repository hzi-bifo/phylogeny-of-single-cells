rule map_reads:
    input:
        reads=expand("results/merged/{{sample}}.{read}.fastq.gz", read=["1", "2"]),
        idx=rules.bwa_index.output,
    output:
        "results/mapped/{sample}.sorted.bam",
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
    resources:
        runtime=179,
        mem_mb=16000,
    threads: 8
    wrapper:
        "v1.21.1/bio/bwa/mem"


rule mark_duplicates:
    input:
        bams="results/mapped/{sample}.sorted.bam",
    output:
        bam=temp("results/markdup/{sample}.sorted.bam"),
        metrics="results/qc/markdup/{sample}.metrics.txt",
    log:
        "logs/picard/markdup/{sample}.log",
    conda:
        "../envs/picard.yaml"
    params:
        java_opts="-Xmx3072m",
        extra="{c} {d}".format(
            c=config["picard"]["markduplicates"],
            d="--TAG_DUPLICATE_SET_MEMBERS true --SORTING_COLLECTION_SIZE_RATIO 0.1",
        ),
        bams=lambda wc, input: list(map("--INPUT {}".format, [ input.bams ] if isinstance(input.bams, str) else input.bams)),
    resources:
        mem_mb=4096,
        runtime=119,
    # TODO: return to using wrapper, once this PR is propagated to the wrapper:
    #   https://github.com/snakemake/snakemake-wrapper-utils/pull/24
    script:
        "../scripts/picard_markduplicates.py"



rule calc_consensus_reads:
    input:
        "results/markdup/{sample}.sorted.bam",
    output:
        consensus_r1=temp("results/consensus/fastq/{sample}.1.fq"),
        consensus_r2=temp("results/consensus/fastq/{sample}.2.fq"),
        consensus_se=temp("results/consensus/fastq/{sample}.se.fq"),
        skipped=temp("results/consensus/{sample}.skipped.bam"),
    log:
        "logs/consensus/{sample}.log",
    params:
        extra="",
    resources:
        mem_mb=lambda wildcards, input: 1.5*input.size_mb,
        runtime=359,
    wrapper:
        "v1.21.1/bio/rbt/collapse_reads_to_fragments-bam"


rule map_consensus_reads:
    input:
        reads=get_processed_consensus_input,
        idx=rules.bwa_index.output,
    output:
        temp("results/consensus/{sample}.consensus.{read_type}.mapped.bam"),
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=lambda w: "-C {}".format(get_read_group(w)),
        sort="samtools",
        sort_order="coordinate",
    resources:
        runtime=29,
        mem_mb=16000,
    wildcard_constraints:
        read_type="pe|se",
    log:
        "logs/bwa_mem/{sample}.{read_type}.consensus.log",
    threads: 8
    wrapper:
        "v1.21.1/bio/bwa/mem"


rule merge_consensus_reads:
    input:
        "results/consensus/{sample}.skipped.bam",
        "results/consensus/{sample}.consensus.se.mapped.bam",
        "results/consensus/{sample}.consensus.pe.mapped.bam",
    output:
        temp("results/consensus/{sample}.merged.bam"),
    log:
        "logs/samtools_merge/{sample}.log",
    threads: 8
    wrapper:
        "v1.21.1/bio/samtools/merge"


rule sort_consensus_reads:
    input:
        "results/consensus/{sample}.merged.bam",
    output:
        temp("results/consensus/{sample}.sorted.bam"),
    log:
        "logs/samtools_sort/{sample}.log",
    threads: 8
    resources:
        mem_mb=lambda wildcards, threads: threads * 900
    wrapper:
        "v1.21.1/bio/samtools/sort"


rule bam_index_consensus:
    input:
        "results/consensus/{sample}.sorted.bam",
    output:
        temp("results/consensus/{sample}.sorted.bai"),
    log:
        "logs/bam_index/consensus/{sample}.sorted.log",
    wrapper:
        "v1.21.1/bio/samtools/index"


rule recalibrate_base_qualities:
    input:
        bam="results/consensus/{sample}.sorted.bam",
        bai="results/consensus/{sample}.sorted.bai",
        ref=rules.create_full_reference.output,
        ref_dict=rules.full_reference_dict.output,
        ref_fai=rules.full_reference_faidx.output,
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table=temp("results/recal/{sample}.grp"),
    conda:
        "../envs/gatk.yaml"
    params:
        extra=config["gatk"]["baserecalibrator"],
        java_opts="-Xmx3072m",
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    threads: 8
    resources:
        mem_mb=4096,
        runtime=59,
    # TODO: return to using wrapper, once this PR is propagated to the wrapper:
    #   https://github.com/snakemake/snakemake-wrapper-utils/pull/24
    script:
        "../scripts/gatk_baserecalibratorspark.py"


rule apply_bqsr:
    input:
        bam="results/consensus/{sample}.sorted.bam",
        bai="results/consensus/{sample}.sorted.bai",
        ref=rules.create_full_reference.output,
        ref_dict=rules.full_reference_dict.output,
        ref_fai=rules.full_reference_faidx.output,
        recal_table="results/recal/{sample}.grp",
    output:
        bam=protected("results/recal/{sample}.sorted.bam"),
        bai="results/recal/{sample}.sorted.bai",
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log",
    conda:
        "../envs/gatk.yaml"
    params:
        extra=config["gatk"]["applybqsr"],  # optional
        java_opts="-Xmx3072m",
    resources:
        mem_mb=4096,
        runtime=59,
    # TODO: return to using wrapper, once this PR is propagated to the wrapper:
    #   https://github.com/snakemake/snakemake-wrapper-utils/pull/24
    script:
        "../scripts/gatk_applybqsr.py"


rule merge_bulks:
    input:
        get_individual_bulk_samples_bam,
    output:
        "results/recal/{individual}.merged_bulk.sorted.bam",
    log:
        "logs/samtools/{individual}.merged_bulk.log",
    params:
        extra="",  # optional additional parameters as string
    threads: 8
    wrapper:
        "v1.21.1/bio/samtools/merge"


rule bam_index_merged_bulks:
    input:
        "results/recal/{individual}.merged_bulk.sorted.bam",
    output:
        "results/recal/{individual}.merged_bulk.sorted.bai",
    log:
        "logs/bam_index/merged_bulks/{individual}.merged_bulk.sorted.log",
    wrapper:
        "v1.21.1/bio/samtools/index"
