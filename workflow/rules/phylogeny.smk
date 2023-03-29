NTS = ["A", "C", "G", "T"]

rule prosolo_probs_to_raxml_ng_ml_gt_and_likelihoods_per_cell:
    """
    This rule can fail with missing output if only very few genomic sites are
    queried, namely if any of the non-homozygous genotypes does not occur in
    the set of called variants. As this is rather unlikelikely to occur, I
    would rather have the rule fail on such cases then to have it produce
    aritificially created empty genotype files.
    """
    input:
        calls="results/calls/{individual}/{sc}.merged_bulk.prosolo.sorted.bcf",
        genotype_mapping=workflow.source_path("../resources/raxml_ng_genotype_mapping.tsv"),
        genotype_order=workflow.source_path("../resources/raxml_ng_genotype_order.csv"),
    output:
        ml=expand(
            "results/raxml_ng_input/{{individual}}/ml_gt_and_likelihoods/{{sc}}_{ref_alt}.tsv",
            ref_alt=[ "_".join([ref, alt]) for ref in NTS for alt in NTS if ref != alt ],
        ),
    log:
        "logs/raxml_ng_input/{individual}/ml_gt_and_likelihoods/{sc}.ref_alt_tsvs.log",
    conda:
        "../envs/vembrane_vlr_bcftools_mlr.yaml"
    params:
        likelihoods_init=lambda wc, input: "$" + "=0.0; $".join(next(csv.reader(open(input.genotype_order)))) + "=0.0;",
        likelihoods_join=lambda wc, input: "$" + ",$".join(next(csv.reader(open(input.genotype_order)))),
        prefix=lambda wc, output: path.dirname(output.ml[0]) + f"/{wc.sc}",
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt * 30 - 1,
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
        "                $MAX = max($HOM_REF, $HET, $HOM_ALT); "
        "                if ($HOM_REF == $MAX) {{ $ONE = $REF; $TWO = $REF; }} "
        "                  elif ($HET == $MAX) {{ $ONE = $REF; $TWO = $ALT }} "
        "                  elif ($HOM_ALT ==$MAX) {{ $ONE = $ALT; $TWO = $ALT }};' "
        "      then join -j REF,ALT -r ONE,TWO --lp ml_ -f {input.genotype_mapping} "
        "      then cut -f CHROM,POS,REF,ALT,ml_IUPAC,likelihoods_{wildcards.sc} "
        "      then rename ml_IUPAC,ml_genotype_{wildcards.sc} "
        "      then reorder -f CHROM,POS,REF,ALT,ml_genotype_{wildcards.sc},likelihoods_{wildcards.sc} "
        "      then split --prefix {params.prefix} -g REF,ALT"
        ") 2> {log}"


# xsv join is more memory efficient, but do it one join at a time
rule join_one_more_cell:
    input:
        sc="results/raxml_ng_input/{individual}/ml_gt_and_likelihoods/{sc}_{ref_alt}.tsv",
        previous_cells="results/raxml_ng_input/{individual}/ml_gt_and_likelihoods/{previous_cells}_{ref_alt}.tsv",
    output:
        "results/raxml_ng_input/{individual}/ml_gt_and_likelihoods/{previous_cells}.{sc}_{ref_alt}.tsv",
    log:
        "logs/raxml_ng_input/{individual}/ml_gt_and_likelihoods/{previous_cells}.{sc}_{ref_alt}.tsv.log",
    conda:
        "../envs/xsv_miller.yaml"
    resources:
        runtime=lambda wildcards, attempt, input: attempt * 30 - 1,
        mem_mb=lambda wildcards, input: input.size_mb * 4,
    threads: 4
    shell:
        "( xsv join --delimiter '\\t' --full CHROM,POS {input.sc} CHROM,POS {input.previous_cells} | " 
        "    xsv fmt --out-delimiter '\\t' | "
        "    mlr --tsv put 'if ( is_empty($CHROM) ) {{ $CHROM = $CHROM_2; $POS = $POS_2 }};' "
        "    then cut -x -f CHROM_2,POS_2,REF_2,ALT_2 "
        "  >{output} "
        ") 2>{log}"

 
rule parse_to_raxml_ng_gt_and_likelihoods:
    input:
        all_cells=lambda wc: expand(
            "results/raxml_ng_input/{{individual}}/ml_gt_and_likelihoods/{cells}_{{ref_alt}}.tsv",
            cells=".".join( samples.loc[ (samples["individual"] == wc.individual) & (samples["sample_type"] == "single_cell"), "sample_name"] )
        ),
    output:
        "results/raxml_ng_input/{individual}/ml_gt_and_likelihoods.{ref_alt}.catg",
    log:
        "logs/raxml_ng_input/{individual}/ml_gt_and_likelihoods.{ref_alt}.catg.log",
    conda:
        "../envs/miller.yaml"
    params:
        default_likelihoods=",".join(["0.1"] * 10),
    threads: 4
    shell:
        "( mlr --from {input.all_cells} --tsv cut -x CHROM,POS "
        "    then unsparsify --fill-with 'N' "
        "    then put "
        '      \'$gt = gsub( joinv( get_values( select($*, func(k,v) {{return k =~ "ml_genotype_.*"}}) ), "," ), ",", "" ); '
        '      for (field, value in select($*, func(k,v) {{return k =~ "likelihoods_.*"}}) ) '
        '        {{ $[field] = ssub(value, "N", "{params.default_likelihoods}" ) }}; '
        '        $without_n = gssub($gt, "N", "");\' '
        '    then filter \'$without_n != ""; '
        '        $without_1 = gsub($without_n, $without_n[1], ""); '
        '        $without_1 != "";\' '
        "    then cut -r -f '^gt$,^likelihoods_.*$' "
        "    then reorder -f gt "
        "  >{output} "
        ") 2>{log}"


rule concat_raxml_ng_input_sites:
    input:
        expand(
            "results/raxml_ng_input/{{individual}}/ml_gt_and_likelihoods.{ref_alt}.catg",
            ref_alt=[ "_".join([ref, alt]) for ref in NTS for alt in NTS if ref != alt ],
        ),
    output:
        "results/raxml_ng_input/{individual}.ml_gt_and_likelihoods.catg",
    log:
        "logs/raxml_ng_input/{individual}.ml_gt_and_likelihoods.catg.log",
    conda:
        "../envs/xsv_sed.yaml"
    threads: 2
    shell:
        "( xsv cat rows --delimiter '\\t' {input} | "
        "    xsv fmt --out-delimiter '\\t' | "
        "    sed -e '1s/gt\\t//' -e '1s/likelihoods_//g' "
        "  >{output}; "
        "  CELLS=$(head -n 1 {output} | wc -w); "
        "  SITES=$(tail -n +2 {output} | wc -l); "
        '  sed -i -e "1i\\$CELLS\\t$SITES" {output} '
        ") 2>{log}"


rule raxml_ng_parse:
    input:
        msa="results/raxml_ng_input/{individual}.ml_gt_and_likelihoods.catg",
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
