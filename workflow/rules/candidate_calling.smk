rule freebayes:
    input:
        ref="resources/reference/full_reference.fa",
        ref_idx="resources/reference/full_reference.fa.fai",
        regions="results/regions/{individual}.covered_regions.filtered.bed",
        # you can have a list of samples here
        samples=lambda w:
            expand("results/recal/{sample}.sorted.bam",
                sample=get_individual_samples(w.individual)
            ),
        index=lambda w:
            expand("results/recal/{sample}.sorted.bai",
                sample=get_individual_samples(w.individual)
            ),

    output:
        "results/candidate-calls/{individual}.freebayes.bcf",
    log:
        "logs/freebayes/{individual}.log",
    params:
        # genotyping is performed by prosolo, hence we deactivate it in freebayes by 
        # always setting --pooled-continuous
        extra="--pooled-continuous --min-alternate-count 1 --min-alternate-fraction {}".format(
            config["freebayes"].get("min_alternate_fraction", "0.01"),
        ),
    threads: 48
    wrapper:
        "v1.2.0/bio/freebayes"
