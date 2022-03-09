rule download_univec:
    output:
        fasta="resources/reference/UniVec{core}/UniVec{core}.fa",
        uv="resources/reference/UniVec{core}/README.uv",
        origins="resources/reference/UniVec{core}/README.vector.origins",
    log:
        "logs/UniVec{core}/download.log",
    conda:
        "../envs/curl.yaml"
    shell:
        "( curl -sS -o {output.fasta} https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec{wildcards.core} && "
        "  curl -sS -o {output.uv} https://ftp.ncbi.nlm.nih.gov/pub/UniVec/README.uv && "
        "  curl -sS -o {output.origins} https://ftp.ncbi.nlm.nih.gov/pub/UniVec/README.vector.origins  ) > {log} 2>&1"


rule download_genome:
    output:
        "resources/reference/genome.fa.gz",
    log:
        "logs/download_genome.log",
    conda:
        "../envs/curl.yaml"
    shell:
        "( curl -sS -o {output} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz "
        "  ) > {log} 2>&1"


rule create_full_reference:
    input:
        genome=rules.download_genome.output,
        univec=get_univec_reference_input,
    output:
        "resources/reference/full_reference.fa"
    log:
        "logs/create_full_reference.log"
    shell:
        "( zcat {input.genome} > {output};"
        "  cat {input.univec} >> {output} "
        " ) > {log} 2>&1"


rule bwa_index:
    input:
        "resources/reference/full_reference.fa",
    output:
        idx=multiext(
            "resources/reference/full_reference",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "logs/bwa_index.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.2.0/bio/bwa/index"
