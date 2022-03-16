rule bam_index:
    input:
        "results/{step}/{sample}.sorted.bam",
    output:
        "results/{step}/{sample}.sorted.bai",
    log:
        "logs/bam_index/{step}/{sample}.sorted.log",
    wrapper:
        "v1.2.0/bio/samtools/index"
