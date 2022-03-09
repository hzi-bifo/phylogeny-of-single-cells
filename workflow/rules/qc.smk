rule fastqc:
    input:
        "{path_and_file_name}.fastq.gz",
    output:
        html="results/qc/fastqc/{path_and_file_name}.html",
        zip="results/qc/fastqc/{path_and_file_name}_fastqc.zip",  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        "--quiet",
    log:
        "logs/fastqc/{path_and_file_name}.log",
    threads: 1
    wrapper:
        "v1.2.0/bio/fastqc"


rule multiqc:
    input:
        get_multiqc_input,
    output:
        "results/qc/multiqc.html",
    params:
        "",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "v1.2.0/bio/multiqc"
