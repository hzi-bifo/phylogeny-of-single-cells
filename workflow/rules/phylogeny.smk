rule prosolo_probs_to_raxml_ng_genotypes_per_cell:
    input:
        "results/genotype_fdr/{individual}/{sc}.merged_bulk.prosolo.sorted.{genotype}.fdr_controlled.bcf",
    output:
        "results/raxml_ng_input/{individual}/per_genotype/{sc}.{genotype}.genotype_likelihoods.tsv",
    log:
        "logs/raxml_ng_input/{individual}/per_genotype/{sc}.{genotype}.genotype_likelihoods.log",
    conda:
        "../envs/vembrane_vlr_bcftools.yaml"
    shell:
        "( varlociraptor decode-phred < {input} |"
        "  vembrane table "
        "    --header 'CHROM, POS, REF, ALT, HOM_REF, HET, HOM_ALT' "
        '    \'CHROM, POS, REF, ALT, INFO["PROB_HOM_REF"] + INFO["PROB_ERR_ALT"], '
        '     INFO["PROB_ADO_TO_ALT"] + INFO["PROB_HET"] + INFO["PROB_ADO_TO_REF"], '
        '     INFO["PROB_HOM_ALT"] + INFO["PROB_ERR_REF"]\' '
        "  >{output}"
        ") 2> {log}"


rule merge_raxml_ng_genotypes_per_cell:
    input:
        hom_ref="results/raxml_ng_input/{individual}/per_genotype/{sc}.hom_ref.genotype_likelihoods.tsv",
        het="results/raxml_ng_input/{individual}/per_genotype/{sc}.het.genotype_likelihoods.tsv",
        hom_alt="results/raxml_ng_input/{individual}/per_genotype/{sc}.hom_alt.genotype_likelihoods.tsv",
    output:
        "results/raxml_ng_input/{individual}/{sc}.genotype_likelihoods.tsv",
    log:
        "logs/raxml_ng_input/{individual}/{sc}.genotype_likelihoods.log",
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/merge_raxml_ng_genotypes_per_cell.R"


rule merge_raxml_ng_genotypes_per_individual:
    input:
        get_all_raxml_gts_for_individual,
    output:
        "results/raxml_ng_input/{individual}.genotype_likelihoods.tsv",
    log:
        "logs/raxml_ng_input/{individual}.genotype_likelihoods.log",
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/merge_raxml_ng_genotypes_per_individual.R"


rule raxml_ng_parse:
    input:
        msa="results/raxml_ng_input/{individual}.genotype_likelihoods.tsv",
    output:
        rba="results/raxml_ng_parse/{individual}.raxml.rba",
        log="logs/raxml_ng_parse/{individual}.raxml.log",
    log:
        "logs/raxml_ng_parse/{individual}.raxml.error.log",
    conda:
        "../envs/raxml_ng.yaml"
    params:
        model=config["raxml_ng"].get("model", "GTGTR+FO"),
        prefix=get_raxml_ng_prefix,
    shell:
        "raxml-ng --parse --msa {input.msa} --model {params.model} --prefix {params.prefix} 2>{log}"


rule raxml_ng:
    input:
        rba="results/raxml_ng_parse/{individual}.raxml.rba",
        log="logs/raxml_ng_parse/{individual}.raxml.log",
    output:
        best_tree="results/raxml_ng/{individual}.raxml.bestTree",
        best_tree_collapsed="results/raxml_ng/{individual}.raxml.bestTreeCollapsed",
        best_model="results/raxml_ng/{individual}.raxml.bestModel",
        bootstraps="results/raxml_ng/{individual}.raxml.bootstraps",
        start_trees="results/raxml_ng/{individual}.raxml.startTree",
        ml_trees="results/raxml_ng/{individual}.raxml.mlTrees",
    log:
        "logs/raxml_ng/{individual}.raxml.error.log",
    conda:
        "../envs/raxml_ng.yaml"
    params:
        model=config["raxml_ng"].get("model", "GTGTR+FO"),
        prefix=get_raxml_ng_prefix,
    threads: get_raxml_ng_threads
    resources:
        mem_mb=get_raxml_ng_mem_mb,
    shell:
        "raxml-ng --all --msa {input.rsa} --model {params.model} --prefix {params.prefix} --threads {threads} --tree pars{{50}},rand{{50}} 2>{log}"


rule raxml_ng_support:
    input:
        best_tree="results/raxml_ng/{individual}.raxml.bestTree",
        bootstraps="results/raxml_ng/{individual}.raxml.bootstraps",
    output:
        support="results/raxml_ng_support/{individual}.raxml.support",
    log:
        "logs/raxml_ng_support/{individual}.raxml.error.log",
    conda:
        "../envs/raxml_ng.yaml"
    params:
        prefix=get_raxml_ng_prefix,
    shell:
        "raxml-ng --support --tree {input.best_tree} --bs-trees {input.bootstraps} --prefix {params.prefix} 2>{log}"


rule raxml_ng_ancestral:
    input:
        rba="results/raxml_ng_parse/{individual}.raxml.rba",
        best_tree="results/raxml_ng/{individual}.raxml.bestTree",
        log="logs/raxml_ng_parse/{individual}.raxml.log",
    output:
        ancestral_tree="results/raxml_ng_ancestral/{individual}.raxml.ancestralTree",
        ancestral_states="results/raxml_ng_ancestral/{individual}.raxml.ancestralStates",
        ancestral_probs="results/raxml_ng_ancestral/{individual}.raxml.ancestralProbs",
    log:
        "logs/raxml_ng_ancestral/{individual}.raxml.error.log",
    conda:
        "../envs/raxml_ng.yaml"
    params:
        model=config["raxml_ng"].get("model", "GTGTR+FO"),
        prefix=get_raxml_ng_prefix,
    threads: get_raxml_ng_threads
    resources:
        mem_mb=get_raxml_ng_mem_mb,
    shell:
        "raxml-ng --ancestral --msa {input.rba} --tree {input.best_tree} --model {params.model} --prefix {params.prefix} --threads {threads} 2>{log}"


rule prosolo_probs_to_scelestial_genotypes_per_cell:
    input:
        "results/genotype_fdr/{individual}/{sc}.merged_bulk.prosolo.sorted.{genotype}.fdr_controlled.bcf",
    output:
        "results/scelestial/{individual}/per_genotype/{sc}.{genotype}.genotypes.tsv",
    log:
        "logs/scelestial/{individual}/per_genotype/{sc}.{genotype}.genotypes.log",
    conda:
        "../envs/vembrane_vlr_bcftools.yaml"
    shell:
        "vembrane table "
        "  --header 'position, REF, ALT "
        "  'f\"{CHROM}:{POS}\", REF, ALT' "
        ">{output} "
        "2> {log} "


rule merge_scelestial_genotypes_per_cell:
    input:
        hom_ref="results/scelestial/{individual}/per_genotype/{sc}.hom_ref.genotypes.tsv",
        het="results/scelestial/{individual}/per_genotype/{sc}.het.genotypes.tsv",
        hom_alt="results/scelestial/{individual}/per_genotype/{sc}.hom_alt.genotypes.tsv",
    output:
        "results/scelestial/{individual}/{sc}.all_genotypes.tsv",
    log:
        "logs/scelestial/{individual}/{sc}.all_genotypes.log",
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/merge_scelestial_genotypes_per_cell.R"


rule merge_scelestial_genotypes_per_individual:
    input:
        get_all_scelestial_gts_for_individual,
    output:
        "results/scelestial/{individual}.genotypes.txt",
    log:
        "logs/scelestial/{individual}.genotypes.log",
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/merge_scelestial_genotypes_per_individual.R"


rule scelestial:
    input:
        "results/scelestial/{individual}.genotypes.txt",
    output:
        "results/scelestial/{individual}.tree.txt",
    log:
        "logs/scelestial/{individual}.tree.log",
    conda:
        "../envs/scelestial.yaml"
    shell:
        "scelestial <{input} >{output} 2>{log}"
