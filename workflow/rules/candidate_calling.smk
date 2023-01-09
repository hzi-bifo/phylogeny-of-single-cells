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
        r" sed 's/\t\([0-9]*\)\t\([0-9]*\)$/:\1-\2/' > {output}"
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
            config["freebayes"].get("min_alternate_fraction", "0.005"),
        ),
    threads: 1
    resources:
        runtime=lambda wildcards, attempt: (12 * attempt - 1) * 60 + 59,
        mem_mb=4990
    shell:
        "(freebayes {params.extra} -r {wildcards.region} -f {input.ref} {input.samples} | "
        " bcftools sort -O b -o {output} -T `mktemp -d` - ) 2> {log}"

rule bcftools_norm_candidate_calls:
    input:
        "results/candidate-calls/{individual}.{region}.freebayes.bcf",
    output:
        "results/candidate-calls/{individual}.{region}.freebayes.norm.bcf",
    log:
        "logs/candidate-calls/{individual}.{region}.freebayes.norm.bcf",
    conda:
        "../envs/bcftools.yaml"
    shell:
        # TODO: turn off the following atomize and instead activate --do-not-normalize
        # once the ProSolo model ist re-integrated with varlociraptor
        "( bcftools norm --multiallelics -any --atomize {input} |"
        # TODO: turn off the following filter to snps once ProSolo is integrated with
        # varlociraptor
        "  bcftools view -O b -o {output} --types snps"
        ") 2>{log}"
