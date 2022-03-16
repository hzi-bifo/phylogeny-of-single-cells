rule falco:
    input:
        "{path_and_file_name}.fastq.gz",
    output:
        "results/qc/falco/{path_and_file_name}/fastqc_data.txt",
        report(
            "results/qc/falco/{path_and_file_name}/fastqc_report.html",
            caption="../report/falco.rst",
            category="Quality Control",
        ),
        "results/qc/falco/{path_and_file_name}/summary.txt",
    params:
        out_dir=lambda wc, output: path.dirname(output[0]),
    conda:
        "../envs/falco.yaml"
    log:
        "logs/falco/{path_and_file_name}.log",
    shell:
        """
        ( falco --outdir {params.out_dir} --dir {params.out_dir} {input} ) >{log} 2>&1
        """


rule multiqc:
    input:
        get_multiqc_input,
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality Control",
        ),
    params:
        "",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "v1.2.0/bio/multiqc"
