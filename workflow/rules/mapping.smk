rule map_reads:
    input:
        reads=expand("results/merged/{{sample}}.{read}.fastq.gz", read=["1", "2"]),
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/{sample}.sorted.bam"),
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "v1.2.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "results/mapped/{sample}.sorted.bam",
    output:
        bam=temp("results/markdup/{sample}.sorted.bam"),
        metrics="results/qc/markdup/{sample}.metrics.txt",
    log:
        "logs/picard/markdup/{sample}.log",
    params:
        extra="{c} {d}".format(
            c=config["picard"]["markduplicates"],
            d="TAG_DUPLICATE_SET_MEMBERS=true",
        ),
        java_opts="-Dpicard.useLegacyParser=false",
    wrapper:
        "v1.2.0/bio/picard/markduplicates"


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
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt collapse-reads-to-fragments bam {input} {output} &> {log}"


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
    wildcard_constraints:
        read_type="pe|se",
    log:
        "logs/bwa_mem/{sample}.{read_type}.consensus.log",
    threads: 8
    wrapper:
        "v1.2.0/bio/bwa/mem"


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
        "v1.2.0/bio/samtools/merge"


rule sort_consensus_reads:
    input:
        "results/consensus/{sample}.merged.bam",
    output:
        temp("results/consensus/{sample}.sorted.bam"),
    log:
        "logs/samtools_sort/{sample}.log",
    threads: 8
    wrapper:
        "v1.2.0/bio/samtools/sort"


rule bam_index:
    input:
        "results/consensus/{sample}.sorted.bam",
    output:
        temp("results/consensus/{sample}.sorted.bai"),
    log:
        "logs/bam_index/consensus/{sample}.sorted.log",
    wrapper:
        "v1.2.0/bio/samtools/index"


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
    params:
        extra=config["gatk"]["baserecalibrator"],
        java_opts="",
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    threads: 8
    wrapper:
        "v1.2.0/bio/gatk/baserecalibratorspark"


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
    params:
        extra=config["gatk"]["applybqsr"],  # optional
        java_opts="",  # optional
    wrapper:
        "v1.2.0/bio/gatk/applybqsr"
