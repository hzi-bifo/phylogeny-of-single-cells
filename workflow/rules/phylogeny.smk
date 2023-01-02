rule prosolo_probs_to_raxml_ng_genotypes_per_cell:
    input:
        "results/genotype_fdr/{individual}/{sc}.merged_bulk.prosolo.sorted.{genotype}.fdr_controlled.bcf",
    output:
        "results/raxml_ng/{individual}/per_genotype/{sc}.{genotype}.genotype_likelihoods.tsv"
    log:
        "logs/raxml_ng/{individual}/per_genotype/{sc}.{genotype}.genotype_likelihoods.log"
    conda:
        "../envs/vembrane_vlr_bcftools.yaml"
    shell:
        "( varlociraptor decode-phred < {input} |"
        "  vembrane table "
        "    --header 'CHROM, POS, REF, ALT, HOM_REF, HET, HOM_ALT' "
        "    'CHROM, POS, REF, ALT, INFO[\"PROB_HOM_REF\"] + INFO[\"PROB_ERR_ALT\"], "
        "     INFO[\"PROB_ADO_TO_ALT\"] + INFO[\"PROB_HET\"] + INFO[\"PROB_ADO_TO_REF\"], "
        "     INFO[\"PROB_HOM_ALT\"] + INFO[\"PROB_ERR_REF\"]' "
        "  >{output}"
        ") 2> {log}"


rule merge_raxml_ng_genotypes_per_cell:
    input:
        hom_ref="results/raxml_ng/{individual}/per_genotype/{sc}.hom_ref.genotype_likelihoods.tsv",
        het="results/raxml_ng/{individual}/per_genotype/{sc}.het.genotype_likelihoods.tsv",
        hom_alt="results/raxml_ng/{individual}/per_genotype/{sc}.hom_alt.genotype_likelihoods.tsv",
    output:
        "results/raxml_ng/{individual}/{sc}.genotype_likelihoods.tsv"
    log:
        "logs/raxml_ng/{individual}/{sc}.genotype_likelihoods.log"
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/merge_raxml_ng_genotypes_per_cell.R"


rule merge_raxml_ng_genotypes_per_individual:
    input:
        get_all_gts_for_individual()
    output:
        "results/raxml_ng/{individual}.genotype_likelihoods.tsv"
    log:
        "logs/raxml_ng/{individual}.genotype_likelihoods.log"
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/merge_raxml_ng_genotypes_per_individual.R"

