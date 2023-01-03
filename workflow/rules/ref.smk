rule download_univec:
    output:
        fasta=temp("resources/reference/UniVec{core}/UniVec{core}.fa"),
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
        temp("resources/reference/genome.fa.gz"),
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
        "resources/reference/full_reference.fa",
    log:
        "logs/create_full_reference.log",
    shell:
        "( zcat {input.genome} | sed -e 's/^>chrM />chrMT /' - | sed -e 's/^>chr/>/' - > {output};"
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
    resources:
        runtime='02:29:00',
        mem_mb=9000,
    wrapper:
        "v1.21.1/bio/bwa/index"


rule full_reference_dict:
    input:
        "resources/reference/full_reference.fa",
    output:
        "resources/reference/full_reference.dict",
    log:
        "logs/samtools/full_reference_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "


rule full_reference_faidx:
    input:
        "resources/reference/full_reference.fa",
    output:
        "resources/reference/full_reference.fa.fai",
    log:
        "logs/samtools/full_reference_faidx.log",
    cache: True
    wrapper:
        "v1.21.1/bio/samtools/faidx"


rule get_known_variants:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="resources/reference/full_reference.fa.fai",
    output:
        vcf="resources/variation.vcf.gz",
    log:
        "logs/ensembl/get_known_variants.log",
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"],
        build=config["ref"]["build"],
        type="all",
    cache: True
    wrapper:
        "v1.21.1/bio/reference/ensembl-variation"


rule remove_iupac_codes:
    input:
        "resources/variation.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz",
    log:
        "logs/rbt/remove_iupac_codes.log",
    resources:
        runtime='00:59:00',
    conda:
        "../envs/rbt.yaml"
    cache: True
    shell:
        "(rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}) 2> {log}"


rule tabix_variation_noiupac:
    input:
        "resources/variation.noiupac.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz.tbi",
    log:
        "logs/tabix/variation.noiupac.vcf.gz.log",
    params:
        "-p vcf",
    cache: True
    wrapper:
        "v1.21.1/bio/tabix/index"
