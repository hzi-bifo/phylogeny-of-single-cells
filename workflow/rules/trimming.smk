rule cutadapt:
    input:
        get_sample_unit_fastqs,
    output:
        fastq1="results/trimmed/{sample}.{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}.{unit}.2.fastq.gz",
        qc="results/trimmed/qc/{sample}.{unit}.qc.txt",
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters=get_cutadapt_parameters,
    resources:
        runtime=lambda wildcards, attempt: 2 * attempt * 60 - 1,
    log:
        "logs/cutadapt/{sample}.{unit}.log",
    threads: 4  # set desired number of threads here
    wrapper:
        "v1.21.1/bio/cutadapt/pe"


rule merge_fastqs:
    input:
        get_merge_fastqs_input,
    output:
        temp("results/merged/{sample}.{read}.fastq.gz"),
    log:
        "logs/merge-fastqs/{sample}.{read}.log",
    resources:
        runtime=39,
    wildcard_constraints:
        read="1|2",
    shell:
        "cat {input} > {output} 2> {log}"
