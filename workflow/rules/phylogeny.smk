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
            "results/raxml_ng/{{individual}}/input/{{individual}}.ml_gt_and_likelihoods/{{sc}}_{ref_alt}.tsv",
            ref_alt=[ "_".join([ref, alt]) for ref in NTS for alt in NTS if ref != alt ],
        ),
    log:
        "logs/raxml_ng/{individual}/input/{individual}.ml_gt_and_likelihoods/{sc}.ref_alt_tsvs.log",
    conda:
        "../envs/vembrane_vlr_bcftools_mlr.yaml"
    params:
        likelihoods_init=lambda wc, input: "$" + "=0.0; $".join(next(csv.reader(open(input.genotype_order)))) + "=0.0;",
        likelihoods_join=lambda wc, input: "$" + ",$".join(next(csv.reader(open(input.genotype_order)))),
        prefix=lambda wc, output: path.dirname(output.ml[0]) + f"/{wc.sc}",
        min_prob_genotype=config.get("min_prob_genotype", 0.98)
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: attempt * 60 * 6 - 1,
    shell:
        # TODO: do bcftools view filtering for {sc} and coverage in in {sc} before prosolo calling, but
        # keep it here for now to avoid rerunning already done prosolo calling
        # background info (TL;DR: we expect linear runtime increases in RAxML-NG
        # with added sites): "Indeed, the amount of work to be done is roughly
        # proportional to the number of patterns for a fixed number of taxa."
        # Pfeiffer and Stamatakis, 2010 ("Hybrid MPI/Pthreads Parallelization of
        # the RAxML Phylogenetics Code")
        "( vembrane filter 'REF != \"N\" and ALT != \"N\"'"
        "    <( bcftools view --samples {wildcards.sc} {input.calls} | "
        "       bcftools view --include 'FORMAT/DP[0]>=4' | "
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
        '                $likelihoods_{wildcards.sc}=joinv([{params.likelihoods_join}], ","); '
        "                $MAX = max($HOM_REF, $HET, $HOM_ALT); "
        '                $clear_evidence_{wildcards.sc} = "N"; '
        "                if ($HOM_REF == $MAX) {{ $ONE = $REF; $TWO = $REF; if ($HOM_REF > {params.min_prob_genotype}) {{ $clear_evidence_{wildcards.sc} = $REF }} }} "
        "                  elif ($HET == $MAX) {{ $ONE = $REF; $TWO = $ALT; if ($HET > {params.min_prob_genotype}) {{ $clear_evidence_{wildcards.sc} = $het_IUPAC }} }} "
        "                  elif ($HOM_ALT ==$MAX) {{ $ONE = $ALT; $TWO = $ALT; if ($HOM_ALT > {params.min_prob_genotype}) {{ $clear_evidence_{wildcards.sc} = $ALT }} }};' "
        "      then join -j REF,ALT -r ONE,TWO --lp ml_ -f {input.genotype_mapping} "
        "      then cut -f CHROM,POS,REF,ALT,clear_evidence_{wildcards.sc},ml_IUPAC,likelihoods_{wildcards.sc} "
        "      then rename ml_IUPAC,ml_genotype_{wildcards.sc} "
        "      then reorder -f CHROM,POS,REF,ALT,clear_evidence_{wildcards.sc},ml_genotype_{wildcards.sc},likelihoods_{wildcards.sc} "
        "      then split --prefix {params.prefix} -g REF,ALT"
        ") 2> {log}"


# xsv join is more memory efficient, but do it one join at a time
rule join_one_more_cell:
    input:
        sc="results/raxml_ng/{individual}/input/{individual}.ml_gt_and_likelihoods/{sc}_{ref_alt}.tsv",
        previous_cells="results/raxml_ng/{individual}/input/{individual}.ml_gt_and_likelihoods/{previous_cells}_{ref_alt}.tsv",
    output:
        temp("results/raxml_ng/{individual}/input/{individual}.ml_gt_and_likelihoods/{previous_cells}.{sc}_{ref_alt}.tsv"),
    log:
        "logs/raxml_ng/{individual}/input/{individual}.ml_gt_and_likelihoods/{previous_cells}.{sc}_{ref_alt}.tsv.log",
    conda:
        "../envs/xsv_miller.yaml"
    resources:
        runtime=lambda wc, input, attempt: attempt * 0.5 * input.size_mb,
        # input.size_mb is only queried for the first input file for the group job, so we
        # need to account for the number of executions of this rule which each adds another
        # single cell to from this individual
        mem_mb=lambda wc, attempt, input: attempt * input.size_mb,
    threads: 2
    shell:
        "( xsv join --delimiter '\\t' --full CHROM,POS {input.sc} CHROM,POS {input.previous_cells} | " 
        "    xsv fmt --out-delimiter '\\t' | "
        "    mlr --tsv put 'if ( is_empty($CHROM) ) {{ $CHROM = $CHROM_2; $POS = $POS_2 }};' "
        "    then cut -x -f CHROM_2,POS_2,REF_2,ALT_2 "
        "    then fill-empty -v 'N' | "
        "  uniq"
        "  >{output} "
        ") 2>{log}"


checkpoint concatenate_cell_mean_coverages:
    input:
        all_cells=lambda wc: expand(
            "results/regions/{{individual}}/{cells}.mosdepth.summary.txt",
            cells=get_single_cells_for_individual(wc.individual)
        ),
    output:
        coverages="results/raxml_ng/input/{individual}.single_cell_coverages.txt"
    log:
        "logs/raxml_ng/input/{individual}.single_cell_coverages.txt"
    conda:
        "../envs/grep.yaml"
    shell:
        'grep "^total" {input.all_cells} >{output.coverages}'


rule parse_to_raxml_ng_gt_and_likelihoods:
    input:
        covered_cells=get_covered_cells_input,
        genotype_mapping=workflow.source_path("../resources/raxml_ng_genotype_mapping.tsv"),
        genotype_order=workflow.source_path("../resources/raxml_ng_genotype_order.csv"),
    output:
        "results/raxml_ng/{individual}/input/{individual}.ml_gt_and_likelihoods.{ref}_{alt}.catg",
    log:
        "logs/raxml_ng/{individual}/input/{individual}.ml_gt_and_likelihoods.{ref}_{alt}.catg.log",
    conda:
        "../envs/miller.yaml"
    params:
        flat_prior=get_flat_prior_for_ref_alt,
    resources:
#        runtime=lambda wc, input, attempt: attempt * 0.5 * input.size_mb,
        mem_mb=20000,
    threads: 8
    shell:
        "( mlr --nr-progress-mod 100000 --from {input.covered_cells} --tsv cut -x -f CHROM,POS,REF,ALT "
        "    then put ' "
        '      $gt = gsub( joinv( get_values( select($*, func(k,v) {{return k =~ "ml_genotype_.*"}}) ), "," ), ",", "" ); '
        '      $clear_evidence = gssub( gsub( joinv( get_values( select($*, func(k,v) {{return k =~ "clear_evidence_.*"}}) ), "," ), ",", "" ), "N", ""); '
        '      for (field, value in select($*, func(k,v) {{return k =~ "likelihoods_.*"}}) ) '
        '        {{ $[field] = ssub(value, "N", "{params.flat_prior}" ) }}; '
        "      ' "
        '    then filter \' $clear_evidence != ""; \' '
        "    then put ' "
        '       $clear_evidence_second = gsub($clear_evidence, $clear_evidence[1], ""); '
        "       ' "
        '    then filter \' $clear_evidence_second != ""; \' '
        "    then cut -r -f '^gt$,^likelihoods_.*$' "
        "    then reorder -f gt "
        "  >{output} "
        ") 2>{log}"


rule filter_to_max_missing_cells:
    input:
        "results/raxml_ng/{individual}/input/{individual}.ml_gt_and_likelihoods.{ref_alt}.catg",
    output:
        temp("results/raxml_ng/{individual}/input/{individual}.max_{n_missing_cells}_missing/ml_gt_and_likelihoods.{ref_alt}.catg"),
    log:
        "logs/raxml_ng/{individual}/input/{individual}.max_{n_missing_cells}_missing/ml_gt_and_likelihoods.{ref_alt}.log",
    conda:
        "../envs/miller.yaml"
    params:
        max_array_length = lambda wc: int(wc.n_missing_cells) + 1,
        mem_mb=20000,
    shell:
        "( mlr --tsv --from {input} filter "
        "   'length(splita($gt, \"N\")) <= {params.max_array_length}'"
        "  >{output}; "
        '  if wc -l {output} | grep -P "^0 "; '
        "  then "
        "    head -n 1 {input} >{output}; "
        "  fi"
        ") 2>{log} "


rule concat_raxml_ng_input_sites:
    input:
        expand(
            "results/raxml_ng/{{individual}}/input/{{individual}}.max_{{n_missing_cells}}_missing/ml_gt_and_likelihoods.{ref_alt}.catg",
            ref_alt=[ "_".join([ref, alt]) for ref in NTS for alt in NTS if ref != alt ],
        ),
    output:
        catg="results/raxml_ng/{individual}/input/max_{n_missing_cells}_missing/{individual}.ml_gt_and_likelihoods.catg",
        variant_sites="results/raxml_ng/{individual}/input/max_{n_missing_cells}_missing/{individual}.ml_gt_and_likelihoods.filtered_sites.txt",
    log:
        "logs/raxml_ng/{individual}/input/max_{n_missing_cells}_missing/{individual}.ml_gt_and_likelihoods.catg.log",
    conda:
        "../envs/xsv_sed.yaml"
    params:
        mem_mb=20000,
    threads: 2
    shell:
        "( xsv cat rows --delimiter '\\t' {input} | "
        "    xsv fmt --out-delimiter '\\t' | "
        "    sed -e '1s/gt\\t//' -e '1s/likelihoods_//g' "
        "  >{output.catg}; "
        "  CELLS=$(head -n 1 {output.catg} | wc -w); "
        "  SITES=$(tail -n +2 {output.catg} | wc -l); "
        "  echo $SITES > {output.variant_sites}; "
        '  sed -i -e "1i\\\\$CELLS\\t$SITES" {output.catg} '
        ") 2>{log}"


rule raxml_ng_parse:
    input:
        msa="results/raxml_ng/{individual}/input/max_{n_missing_cells}_missing/{individual}.ml_gt_and_likelihoods.catg",
    output:
        "results/raxml_ng/{individual}/parse/{model}/max_{n_missing_cells}_missing/{individual}.raxml.log",
    log:
        "logs/raxml_ng/{individual}/parse/{model}/max_{n_missing_cells}_missing/{individual}.raxml.error.log",
    conda:
        "../envs/raxml_ng.yaml"
    params:
        prefix=get_raxml_ng_prefix,
    resources:
        runtime=lambda wildcards, attempt: attempt * 60 * 2 - 1,
        mem_mb=lambda wildcards, attempt, input: attempt * 10 * input.size_mb,
    threads: 2
    shell:
        "(raxml-ng --parse "
        "  --msa {input.msa} "
        "  --model {wildcards.model} "
        "  --prefix {params.prefix} "
        ") 2>{log}"


rule raxml_ng_tree_search:
    input:
        msa="results/raxml_ng/{individual}/input/max_{n_missing_cells}_missing/{individual}.ml_gt_and_likelihoods.catg",
        log="results/raxml_ng/{individual}/parse/{model}/max_{n_missing_cells}_missing/{individual}.raxml.log",
    output:
        best_tree="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.bestTree",
        ml_trees="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.mlTrees", # is not produced with a single starting tree
#        best_tree_collapsed="results/raxml_ng/{individual}.raxml.bestTreeCollapsed", # is not always produced
        best_model="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.bestModel",
        start_trees="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.startTree",
        log="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.log",
    log:
        "logs/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.error.log",
    conda:
        "../envs/raxml_ng.yaml"
    params:
        prefix=get_raxml_ng_prefix,
#    threads: get_raxml_ng_threads
    threads: 32
    resources:
        mem_mb=get_raxml_ng_mem_mb,
        runtime=lambda wildcards, attempt: attempt * 60 * 24 * 2,
    shell:
        "raxml-ng --search "
        "  --msa {input.msa} "
        "  --model {wildcards.model} "
        "  --blmin 1e-9 "
        "  --prefix {params.prefix} "
        "  --prob-msa on "
        "  --threads auto{{{threads}}} "
        "  --tree pars{{250}},rand{{250}} "
        "  --extra pars-par " # compute parsimony trees in parallel
        "  --redo "
        "2>{log}"


rule raxml_ng_bootstrap:
    input:
        msa="results/raxml_ng/{individual}/input/max_{n_missing_cells}_missing/{individual}.ml_gt_and_likelihoods.catg",
        log="results/raxml_ng/{individual}/parse/{model}/max_{n_missing_cells}_missing/{individual}.raxml.log",
    output:
        bootstraps="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bootstraps.raxml.bootstraps",
        log="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bootstraps.raxml.log",
    log:
        "logs/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bootstraps.raxml.error.log",
    conda:
        "../envs/raxml_ng.yaml"
    params:
        prefix=get_raxml_ng_prefix,
#    threads: get_raxml_ng_threads
    threads: 32
    resources:
        mem_mb=get_raxml_ng_mem_mb,
        runtime=lambda wildcards, attempt: attempt * 60 * 24 * 3,
    shell:
        "raxml-ng --bootstrap "
        "  --msa {input.msa} "
        "  --model {wildcards.model} "
        "  --blmin 1e-9 "
        "  --prefix {params.prefix} "
        "  --prob-msa on "
        "  --threads auto{{{threads}}} "
        "  --extra pars-par " # compute parsimony trees in parallel
        "  --bs-trees 2000 "
        "  --extra bs-start-rand " # try avoiding parsimony bias in tree topology bootstrapping
        "  --redo "
        "2>{log}"


rule raxml_ng_bsconverge:
    input:
        bootstraps="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bootstraps.raxml.bootstraps",
        log="results/raxml_ng/{individual}/parse/{model}/max_{n_missing_cells}_missing/{individual}.raxml.log",
    output:
        log="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bsconverge.raxml.log",
    log:
        "logs/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bsconverge.raxml.error.log",
    conda:
        "../envs/raxml_ng.yaml"
    params:
        prefix=get_raxml_ng_prefix,
#    threads: get_raxml_ng_threads
    threads: 4
    resources:
        mem_mb=get_raxml_ng_mem_mb,
        runtime=lambda wildcards, attempt: attempt * 60 * 24 * 3,
    shell:
        "raxml-ng --bsconverge "
        "  --bs-trees {input.bootstraps} "
        "  --prefix {params.prefix} "
        "  --threads auto{{{threads}}} "
        "  --redo "
        "2>{log}"


rule raxml_ng_support:
    input:
        msa="results/raxml_ng/{individual}/input/max_{n_missing_cells}_missing/{individual}.ml_gt_and_likelihoods.catg",
        log="results/raxml_ng/{individual}/parse/{model}/max_{n_missing_cells}_missing/{individual}.raxml.log",
        best_tree="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.bestTree",
        bootstraps="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bootstraps.raxml.bootstraps",
    output:
        bootstraps="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.support.raxml.support",
        log="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.support.raxml.log",
    log:
        "logs/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.support.raxml.error.log",
    conda:
        "../envs/raxml_ng.yaml"
    params:
        prefix=get_raxml_ng_prefix,
#    threads: get_raxml_ng_threads
    threads: 4
    resources:
        mem_mb=get_raxml_ng_mem_mb,
        runtime=lambda wildcards, attempt: attempt * 60 * 24 * 3,
    shell:
        "raxml-ng --support "
        "  --tree {input.best_tree} "
        "  --bs-trees {input.bootstraps} "
        "  --prefix {params.prefix} "
        "  --bs-metric tbe "
        "  --redo "
        "2>{log}"


rule raxml_ng_rfdist_across_trees:
    input:
        trees=lambda wc: expand(
            "results/raxml_ng/{{individual}}/results/{{model}}/max_{{n_missing_cells}}_missing/{{individual}}.{{model}}.max_{{n_missing_cells}}_missing.{infix}.raxml.{{type}}",
            infix="bootstraps" if wc.type == "bootstraps" else "search",
        ),
    output:
        rf_dist="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.rf_dist.raxml.rfDistances",
        log="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.rf_dist.raxml.log",
    log:
        "logs/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.rf_dist.raxml.rfDistances.log",
    conda:
        "../envs/raxml_ng.yaml"
    params:
        prefix=get_raxml_ng_prefix,
    threads: 2
    shell:
        "raxml-ng --rfdist "
        "  --tree {input.trees} "
        "  --prefix {params.prefix} "
        "  --redo "
        "2>{log}"

rule raxml_ng_rfdist_to_tsv:
    input:
        rf_dist="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.raxml.rfDistances",
    output:
        tsv="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.raxml.rf_dist.tsv",
    log:
        "logs/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.raxml.rfDistances_to_tsv.log",
    run:
        with open(input.rf_dist) as f:
            for line in f:
                if line.startswith("Loaded "):
                    total_trees = line.split("Loaded ")[1].split(" trees with ")[0]
                elif line.startswith("Average absolute RF distance in this tree set:"):
                    absolute_rf_distance = line.split("in this tree set: ")[1]
                elif line.startswith("Average relative RF distance in this tree set:"):
                    relative_rf_distance = line.split("in this tree set: ")[1]
                elif line.startswith("Number of unique topologies in this tree set:"):
                    unique_topologies = line.split("in this tree set: ")[1]
        with open(output.tsv) as o:
            o.write("\t".join(
                [
                    "inidividual",
                    "model",
                    "max_missing",
                    "tree_type",
                    "total_trees",
                    "absolute_rf_distance",
                    "relative_rf_distance",
                    "unique_topologies",
                ]
            ))
            o.write("\t".join(
                [
                    wildcards.inidividual,
                    wildcards.model,
                    wildcards.max_missing,
                    wildcards.type,
                    total_trees,
                    absolute_rf_distance,
                    relative_rf_distance,
                    unique_topologies,
                ]
            ))


rule concatenate_rf_dist_results_per_individual:
    input:
        rf_dists=expand(
            "results/raxml_ng/{{individual}}/results/{model}/max_{n_missing_cells}_missing/{{individual}}.{model}.max_{n_missing_cells}_missing.{type}.raxml.rf_dist.tsv",
            model=config["raxml_ng"]["models"],
            n_missing_cells=config["raxml_ng"]["max_missing"],
            type=["startTree", "mlTrees", "bootstraps"],
        ),
    output:
        tsv="results/trees/{individual}.max_missing_stable_topology_selection.tsv",
    log:
        "logs/trees/{individual}.max_missing_stable_topology_selection.log",
    conda:
        "../envs/xsv.yaml"
    shell:
        "( xsv cat rows "
        '    --delimiter "\\t" '
        "    {input.rf_dists} | "
        "  xsv fmt --out-delimiter '\\t' "
        "  >{output.tsv} "
        ") 2>{log}"


rule plot_distinct_topologies_across_missingness:
    input:
        tsv="results/trees/{individual}.max_missing_stable_topology_selection.tsv",
    output:
        plot="results/trees/{individual}.max_missing_stable_topology_selection.pdf",
    log:
        "logs/trees/{individual}.max_missing_stable_topology_selection.log"
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/plot_max_missing_topologies.R"


rule plot_support_tree:
    input:
        support="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.support.raxml.support",
        samples=config["samples"],
    output:
        support="results/trees/{individual}/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.raxml.support.pdf",
    log:
        "logs/trees/{model}/max_{n_missing_cells}_missing/{individual}.raxml.support.log",
    conda:
        "../envs/ggtree.yaml"
    script:
        "../scripts/plot_support_tree.R"


rule raxml_ng_ancestral:
    input:
        msa="results/raxml_ng/input/{individual}.ml_gt_and_likelihoods.catg",
        best_tree="results/raxml_ng/{individual}.raxml.bestTree",
        log="results/raxml_ng/parse/{individual}.raxml.log",
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
        "(raxml-ng --ancestral "
        "  --msa {input.msa} "
        "  --tree {input.best_tree} "
        "  --model {params.model} "
        "  --prefix {params.prefix} "
        "  --threads {threads} "
        ") 2>{log}"


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
