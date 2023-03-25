rule prosolo_probs_to_raxml_ng_ml_gt_and_likelihoods_per_cell:
    input:
        calls="results/calls/{individual}/{sc}.merged_bulk.prosolo.sorted.bcf",
        genotype_mapping=workflow.source_path("../resources/raxml_ng_genotype_mapping.tsv"),
        genotype_order=workflow.source_path("../resources/raxml_ng_genotype_order.csv"),
    output:
        ml="results/raxml_ng_input/{individual}/{sc}.ml_gt_and_likelihoods.tsv",
    log:
        "logs/raxml_ng_input/{individual}/{sc}.ml_gt_and_likelihoods.log",
    conda:
        "../envs/vembrane_vlr_bcftools_mlr.yaml"
    params:
        likelihoods_init=lambda wc, input: "$" + "=0.0; $".join(next(csv.reader(open(input.genotype_order)))) + "=0.0;",
        likelihoods_join=lambda wc, input: "$" + ",$".join(next(csv.reader(open(input.genotype_order)))),
    threads: 6
    resources:
        runtime=lambda wildcards, attempt: attempt * 60 - 1,
    shell:
        # TODO: do bcftools view filtering for {sc} and coverage in in {sc} before prosolo calling, but
        # keep it here for now to avoid rerunning already done prosolo calling
        "( vembrane filter 'REF != \"N\" and ALT != \"N\"'"
        "    <( bcftools view --samples {wildcards.sc} {input.calls} | "
        "       bcftools view --include 'FORMAT/DP[0]>=1' | "
        "       varlociraptor decode-phred "
        "     ) |"
        "  vembrane table "
        "    --header 'CHROM, POS, REF, ALT, HOM_REF, HET, HOM_ALT' "
        '    \'CHROM, POS, REF, ALT, INFO["PROB_HOM_REF"] + INFO["PROB_ERR_ALT"], '
        '     INFO["PROB_ADO_TO_ALT"] + INFO["PROB_HET"] + INFO["PROB_ADO_TO_REF"], '
        '     INFO["PROB_HOM_ALT"] + INFO["PROB_ERR_REF"]\' | '
        "  mlr --tsv join -j REF,ALT --lp het_ -f {input.genotype_mapping} "
        "      then put '{params.likelihoods_init} "
        "                $[$REF] = $HOM_REF; "
        "                $[$het_IUPAC] = $HET; "
        "                $[$ALT] = $HOM_ALT; "
        "                $likelihoods_{wildcards.sc}=joinv([{params.likelihoods_join}], \",\"); "
        "                $variant_key = format(\"{{}}:{{}}_{{}}_{{}}\", $CHROM,$POS,$REF,$ALT); "
        "                $MAX = max($HOM_REF, $HET, $HOM_ALT); "
        "                if ($HOM_REF == $MAX) {{ $ONE = $REF; $TWO = $REF; }} "
        "                  elif ($HET == $MAX) {{ $ONE = $REF; $TWO = $ALT }} "
        "                  elif ($HOM_ALT ==$MAX) {{ $ONE = $ALT; $TWO = $ALT }};' "
        "      then join -j REF,ALT -r ONE,TWO --lp ml_ -f {input.genotype_mapping} "
        "      then cut -f variant_key,ml_IUPAC,likelihoods_{wildcards.sc} "
        "      then rename ml_IUPAC,ml_genotype_{wildcards.sc} "
        "      then reorder -f variant_key,ml_genotype_{wildcards.sc},likelihoods_{wildcards.sc} "
        "  >{output.ml}"
        ") 2> {log}"


rule merge_raxml_ng_genotypes_per_individual:
    input:
        get_all_raxml_gts_for_individual,
    output:
        "results/raxml_ng_input/{individual}.ml_gt_and_likelihoods.tsv",
    log:
        "logs/raxml_ng_input/{individual}.ml_gt_and_likelihoods.log",
    conda:
        "../envs/csvkit.yaml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 90 - 1,
        mem_mb=lambda wildcards, input: input.size_mb * 1.1,
    shell:
        "( csvjoin --outer "
        "    --tabs "
        "    --out-tabs "
        "    --quoting 3 " # 3 = no quoting
        "    --no-inference " # major bottleneck according to: https://csvkit.readthedocs.io/en/1.1.1/tricks.html#slow-performance
        "    --columns variant_key "
        "    {input} "
        "    >{output} "
        ") 2>{log}"


rule raxml_ng_parse:
    input:
        msa="results/raxml_ng_input/{individual}.ml_gt_and_likelihoods.tsv",
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
        "raxml-ng --all --prob-msa --msa {input.rsa} --model {params.model} --prefix {params.prefix} --threads {threads} --tree pars{{50}},rand{{50}} 2>{log}"


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
        "../envs/vembrane_vlr.yaml"
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
