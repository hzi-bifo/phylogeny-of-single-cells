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
    log:
        "logs/cutadapt/{sample}.{unit}.log",
    threads: 4  # set desired number of threads here
    wrapper:
        "v1.1.0/bio/cutadapt/pe"


rule merge_fastqs:
    input:
        expand(
            "results/trimmed/{{sample}}.{unit}.{{read}}.fastq.gz",
            unit=lambda wc: units.loc[wc.sample, "unit_name"]
        )
    output:
        temp("results/merged/{sample}.{read}.fastq.gz"),
    log:
        "logs/merge-fastqs/{sample}.{read}.log",
    wildcard_constraints:
        read="1|2",
    shell:
        "cat {input} > {output} 2> {log}"
