checkpoint create_freebayes_regions:
    input:
        ref="resources/reference/full_reference.fa",
        ref_idx="resources/reference/full_reference.fa.fai",
        regions="results/regions/{individual}.covered_regions.filtered.bed",
    output:
        "results/regions/{individual}.freebayes_regions.bed",
    log:
        "logs/regions/{individual}.freebayes_regions.log",
    conda:
        "../envs/freebayes.yaml"
    params:
        chunksize=lambda w: config["freebayes"]["chunksize"],
    shell:
        "( bedtools intersect -a "
        r"  <(sed 's/:\([0-9]*\)-\([0-9]*\)$/\t\1\t\2/' "
        "     <(fasta_generate_regions.py {input.ref}.fai {params.chunksize} )"
        "    ) "
        "   -b {input.regions} | "
        r" sed 's/\t\([0-9]*\)\t\([0-9]*\)$/\/\1-\2/' > {output}"
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
        "results/candidate_calls/{individual}/{chromosome}/{region}.freebayes.bcf",
    log:
        "logs/freebayes/{individual}/{chromosome}/{region}.log",
    conda:
        "../envs/freebayes.yaml"
    params:
        # genotyping is performed by prosolo, hence we deactivate it in freebayes by 
        # always setting --pooled-continuous
        extra="--pooled-continuous --use-best-n-alleles 4 --min-alternate-count 2 --min-alternate-total 4 --min-alternate-fraction {}".format(
            config["freebayes"].get("min_alternate_fraction", "0.005"),
        ),
    threads: 1
    resources:
        runtime=lambda wildcards, attempt: 5 * attempt * 60 - 1,
        mem_mb=lambda wc, attempt: 8000 + 16000 * (attempt - 1),
    shell:
        "(freebayes {params.extra} -r {wildcards.chromosome}:{wildcards.region} -f {input.ref} {input.samples} | "
        " bcftools sort -O b -o {output} -T `mktemp -d` - ) 2> {log}"


rule bcftools_norm_candidate_calls:
    input:
        "results/candidate_calls/{individual}/{chromosome}/{region}.freebayes.bcf",
    output:
        "results/candidate_calls/{individual}/{chromosome}/{region}.freebayes.norm.bcf",
    log:
        "logs/candidate_calls/{individual}/{chromosome}/{region}.freebayes.norm.bcf",
    conda:
        "../envs/bcftools.yaml"
    resources:
        runtime=lambda wildcards, attempt: 30 * attempt - 1,
    shell:
        # TODO: turn off the following atomize and instead activate --do-not-normalize
        # once the ProSolo model ist re-integrated with varlociraptor
        # TODO: turn off the following filter to snps once ProSolo is integrated with
        # varlociraptor
        "( bcftools norm --multiallelics -any --atomize {input} |"
        "  bcftools view -O b -o {output} --types snps"
        ") 2>{log}"


rule bcftools_index_candidate_calls:
    input:
        "results/candidate_calls/{individual}/{chromosome}/{region}.freebayes.norm.bcf",
    output:
        "results/candidate_calls/{individual}/{chromosome}/{region}.freebayes.norm.bcf.csi",
    log:
        "logs/candidate_calls/bcftools_index/{individual}/{chromosome}/{region}.freebayes.norm.log",
    params:
        extra="",  # optional parameters for bcftools index
    resources:
        runtime=lambda wildcards, attempt: 20 * attempt - 1,
    wrapper:
        "v1.22.0/bio/bcftools/index"


rule aggregate_freebayes_region_calls:
    input:
        calls=aggregate_freebayes_region_calls_input(),
        indexes=aggregate_freebayes_region_calls_input(ext=".bcf.csi"),
    output:
        "results/candidate_calls/{individual}.freebayes.norm.bcf",
    log:
        "logs/candidate_calls/{individual}.freebayes.norm.log",
    params:
        extra="-a",
    resources:
        runtime=lambda wildcards, attempt: 6 * 60 * attempt - 1,
        mem_mb=lambda wildcards, input, attempt: input.size_mb * 4 * attempt
    wrapper:
        "v1.21.1/bio/bcftools/concat"
