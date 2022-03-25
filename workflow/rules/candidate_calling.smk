checkpoint create_freebayes_regions:
    input:
        ref="resources/reference/full_reference.fa",
        ref_idx="resources/reference/full_reference.fa.fai",
        regions="results/regions/{individual}.covered_regions.filtered.bed",
    output:
        "results/regions/{individual}.freebayes_regions.bed"
    log:
        "logs/regions/{individual}.freebayes_regions.log"
    conda:
        "../envs/freebayes.yaml"
    params:
        chunksize=lambda w: config["freebayes"]["chunksize"]
    shell:
        "( bedtools intersect -a "
        r"  <(sed 's/:\([0-9]*\)-\([0-9]*\)$/\t\1\t\2/' "
        "     <(fasta_generate_regions.py {input.ref}.fai {params.chunksize} )"
        "    ) "
        "   -b {input.regions} | "
        r" sed 's/\t\([0-9]*\)\t\([0-9]*\)$/:\1-\2/') > {output}"
        "  ) 2> {log} "


rule freebayes_per_region:
    input:
        ref="resources/reference/full_reference.fa",
        ref_idx="resources/reference/full_reference.fa.fai",
        # you can have a list of samples here
        samples=lambda w: expand(
            "results/recal/{sample}.sorted.bam",
            sample=get_individual_samples(w.individual),
        ),
        index=lambda w: expand(
            "results/recal/{sample}.sorted.bai",
            sample=get_individual_samples(w.individual),
        ),
    output:
        "results/candidate-calls/{individual}.{region}.freebayes.bcf",
    log:
        "logs/freebayes/{individual}.{region}.log",
    conda:
        "../envs/freebayes.yaml"
    params:
        # genotyping is performed by prosolo, hence we deactivate it in freebayes by 
        # always setting --pooled-continuous
        extra="--pooled-continuous --min-alternate-count 1 --min-alternate-total 2 --min-alternate-fraction {}".format(
            config["freebayes"].get("min_alternate_fraction", "0.01"),
        ),
    threads: 2
    shell:
        "(freebayes {params.extra} -f {input.ref} {input.samples} | "
        " bcftools sort -O b -o {output} - ) 2> {log}"


rule aggregate_freebayes:
    input:
        aggregate_freebayes_input
    output:
        "results/candidate-calls/{individual}.freebayes.bcf",
    log:
        "logs/aggregate_candidate_calls/{individual}.freebayes.bcf",
    params:
        uncompressed_bcf=False,
        extra="--allow-overlaps",  # optional parameters for bcftools concat (except -o)
    threads: 4
    resources:
        mem_mb=16000,
    wrapper:
        "v1.3.1/bio/bcftools/concat"


rule scatter_candidates:
    input:
        "results/candidate-calls/{individual}.freebayes.bcf",
    output:
        scatter.calling(
            "results/candidate-calls/{{individual}}.freebayes.{scatteritem}.bcf"
        ),
    log:
        "logs/scatter_candidates/{individual}.freebayes.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"
