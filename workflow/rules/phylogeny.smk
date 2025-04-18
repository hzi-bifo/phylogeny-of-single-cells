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
    threads: 16 # it can make sense to set this to --set-threads raxml_ng_bootstrap=1 for jobs that are failing, because parallel execution sometimes jumbles up the log file order that we need to parse out likelihoods downstream
    resources:
        mem_mb=get_raxml_ng_mem_mb,
        runtime=lambda wildcards, attempt: attempt * 60 * 22,
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


rule gotree_collapse_bootstrap_trees:
    input:
        bootstraps="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bootstraps.raxml.bootstraps",
        variant_sites="results/raxml_ng/{individual}/input/max_{n_missing_cells}_missing/{individual}.ml_gt_and_likelihoods.filtered_sites.txt",
    output:
        collapsed="results/gotree/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bootstraps.raxml.collapsed.bootstraps",
    log:
        "logs/gotree/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bootstraps.raxml.collapsed.bootstraps.log",
    conda:
        "../envs/gotree.yaml"
    params:
        min_length=lambda wc, input: min(1.0e-04, 0.1/pd.read_csv(input.variant_sites, header=None).loc[0,0])
    shell:
        "gotree collapse length --length {params.min_length} --input {input.bootstraps} --output {output.collapsed} 2>{log} "


rule gotree_collapse_tree:
    input:
        best_tree="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.{type}",
        variant_sites="results/raxml_ng/{individual}/input/max_{n_missing_cells}_missing/{individual}.ml_gt_and_likelihoods.filtered_sites.txt",
    output:
        collapsed="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.collapsed.{type}",
    log:
        "logs/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.collapsed.{type}.log",
    conda:
        "../envs/gotree.yaml"
    params:
        min_length=lambda wc, input: min(1.0e-04, 0.1/pd.read_csv(input.variant_sites, header=None).loc[0,0])
    shell:
        "gotree collapse length --length {params.min_length} --input {input.best_tree} --output {output.collapsed} 2>{log} "


rule extract_ml_tree_likelihoods:
    input:
        log="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{infix}.raxml.log",
    output:
        tsv="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{infix}.raxml.log_likelihoods.tsv",
    log:
        "logs/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{infix}.raxml.log_likelihoods.log",
    conda:
        "../envs/grep.yaml"
    params:
        search_string=lambda wc: "Bootstrap tree" if wc.infix == "bootstraps" else "ML tree search"
    shell:
        '( grep "{params.search_string} #" {input.log} | '
        r"sed -r -e 's/^.*{params.search_string} #([0-9]+), logLikelihood: (-[0-9]+.[0-9]+)$/\1\t\2/' | "
        "sort -n -k 1,1 "
        ">{output.tsv} "
        ") 2>{log}"


rule cluster_info_dist_across_trees:
    input:
        trees=lambda wc: expand(
            "results/{{software}}/{{individual}}/results/{{model}}/max_{{n_missing_cells}}_missing/{{individual}}.{{model}}.max_{{n_missing_cells}}_missing.{infix}.raxml.collapsed.{{type}}",
            infix="bootstraps" if wc.type == "bootstraps" else "search",
        ),
        likelihoods=lambda wc: expand(
            "results/raxml_ng/{{individual}}/results/{{model}}/max_{{n_missing_cells}}_missing/{{individual}}.{{model}}.max_{{n_missing_cells}}_missing.{infix}.raxml.log_likelihoods.tsv",
            infix="bootstraps" if wc.type == "bootstraps" else "search",
        )
    output:
        all_cid="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.all.cluster_info_dist.tsv.gz",
        silhouette_scores="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.all.clustering_silhouette_scores.pdf",
        best_clustering_log_likelihoods="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.all.best_clustering_log_likelihoods.pdf",
        best_clustering_silhouette_scores="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.all.best_clustering_silhouette_scores.pdf",
        best_clustering_silhouette_ll="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.all.best_clustering.silhouette_scores_log_likelihoods.tsv",
        best_clustering_stats="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.all.best_clustering.stats.tsv",
        mapping_cluster_likelihood="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.all.best_clustering_pcoa_mapping_with_tree_likelihoods.pdf",
        best_cluster_cid="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.best_cluster.best_cluster_info_dist.tsv",
        best_cluster_trees="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.best_cluster.best_cluster_trees.newick",
    log:
        "logs/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.cluster_info_dist.log",
    conda:
        "../envs/treedist.yaml"
    threads: 1
    script:
        "../scripts/cluster_info_dist.R"


rule gotree_compute_consensus_tree:
    input:
        best_cluster_trees="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.mlTrees.best_cluster.best_cluster_trees.newick",
    output:
        consensus_tree="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.consensusTree",
    log:
        "logs/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.consensusTree.log",
    conda:
        "../envs/gotree.yaml"
    threads: 4
    shell:
        "( gotree compute consensus "
        "    --input {input.best_cluster_trees} "
        "    --freq-min 0.5 " # this should be the default, but without a min frequency we get an error about this
        "    --threads {threads} "
        "    --output {output.consensus_tree} "
        ") 2>{log}"


rule gotree_support_collapsed_trees:
    input:
        ref_tree="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.collapsed.{tree_type}",
        bootstraps="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bootstraps.raxml.collapsed.bootstraps",
    output:
        support="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.support.{tree_type}.raxml.collapsed.support",
    log:
        "logs/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.support.{tree_type}.raxml.support.log",
    conda:
        "../envs/gotree.yaml"
    threads: 4
    shell:
        "( gotree compute support tbe "
        "    --threads {threads} "
        "    --bootstrap {input.bootstraps} "
        "    --reftree {input.ref_tree} "
        "    --out {output.support} "
        ") 2>{log}"


rule raxml_ng_rfdist_across_trees:
    input:
        trees=lambda wc: expand(
            "results/{{software}}/{{individual}}/results/{{model}}/max_{{n_missing_cells}}_missing/{{individual}}.{{model}}.max_{{n_missing_cells}}_missing.{infix}.raxml.{{type}}",
            infix="bootstraps" if wc.type == "bootstraps" else "search",
        ),
    output:
        rf_dist="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.rf_dist.raxml.rfDistances",
        log="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.rf_dist.raxml.log",
    log:
        "logs/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.rf_dist.raxml.rfDistances.log",
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
        rf_dist="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.rf_dist.raxml.log",
    output:
        tsv="results/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.raxml.rf_dist.tsv",
    log:
        "logs/{software}/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.{type}.rf_dist.raxml.rfDistances_to_tsv.log",
    run:
        with open(input.rf_dist) as f:
            for line in f:
                if line.startswith("Loaded "):
                    total_trees = line.split("Loaded ")[1].split(" trees with ")[0]
                elif line.startswith("Average absolute RF distance in this tree set:"):
                    absolute_rf_distance = line.split("in this tree set: ")[1].rstrip()
                elif line.startswith("Average relative RF distance in this tree set:"):
                    relative_rf_distance = line.split("in this tree set: ")[1].rstrip()
                elif line.startswith("Number of unique topologies in this tree set:"):
                    unique_topologies = line.split("in this tree set: ")[1].rstrip()
        with open(output.tsv, mode="w") as o:
            o.write("\t".join(
                [
                    "individual",
                    "model",
                    "max_missing",
                    "tree_type",
                    "total_trees",
                    "absolute_rf_distance",
                    "relative_rf_distance",
                    "unique_topologies",
                ]
            ))
            o.write("\n")
            o.write("\t".join(
                [
                    wildcards.individual,
                    wildcards.model,
                    wildcards.n_missing_cells,
                    wildcards.type,
                    total_trees,
                    absolute_rf_distance,
                    relative_rf_distance,
                    unique_topologies,
                ]
            ))



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
        runtime=lambda wildcards, attempt: attempt * 60 * 22,
    shell:
        "raxml-ng --bsconverge "
        "  --bs-trees {input.bootstraps} "
        "  --prefix {params.prefix} "
        "  --threads auto{{{threads}}} "
        "  --redo "
        "2>{log}"


rule raxml_ng_bsconverge_to_tsv:
    input:
        bsconverge="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bsconverge.raxml.log",
        ml_search="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.search.raxml.log",
    output:
        tsv="results/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bootstraps.raxml.bsconverge.tsv",
    log:
        "logs/raxml_ng/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.bootstraps.bsconverge.raxml.log",
    run:
        with open(input.bsconverge) as f:
            for line in f:
                if line.startswith("Loaded "):
                    total_trees = line.split("Loaded ")[1].split(" trees with ")[0]
                elif line.startswith("Bootstopping test converged after "):
                    convergence = line.split("Bootstopping test converged after ")[1].split(" trees")[0]
                elif line.startswith("Bootstopping test did not converge after "):
                    convergence = total_trees
        with open(input.ml_search) as f:
            for line in f:
                if line.startswith("Alignment sites: "):
                    usable_sites = line.split("Alignment sites: ")[1].rstrip()
                elif line.startswith("AIC score: "):
                    scores = line.rstrip().split(" / ")
                    [ aic, aicc, bic ] = [ s.split(" score: ")[1] for s in scores]
                    break
        with open(output.tsv, mode="w") as o:
            o.write("\t".join(
                [
                    "individual",
                    "model",
                    "max_missing",
                    "usable_sites",
                    "tree_type",
                    "total_trees",
                    "convergence",
                    "aic_score",
                    "aicc_score",
                    "bic_score",
                ]
            ))
            o.write("\n")
            o.write("\t".join(
                [
                    wildcards.individual,
                    wildcards.model,
                    wildcards.n_missing_cells,
                    usable_sites,
                    "bootstraps",
                    total_trees,
                    convergence,
                    aic,
                    aicc,
                    bic,
                ]
            ))


rule concatenate_raxml_ng_tsvs_per_individual_per_metric:
    input:
        rf_dists=lambda wc: expand(
            "results/raxml_ng/{{individual}}/results/{model}/max_{n_missing_cells}_missing/{{individual}}.{model}.max_{n_missing_cells}_missing.{type}.raxml.{{metric}}.tsv",
            model=config["raxml_ng"]["models"],
            n_missing_cells=config["raxml_ng"]["max_missing"],
            type=["startTree", "mlTrees", "bootstraps"] if wc.metric == "rf_dist" else ["bootstraps"],
        ),
    output:
        tsv="results/trees/{individual}.{metric}.max_missing_stable_topology_selection.tsv",
    log:
        "logs/trees/{individual}.{metric}.max_missing_stable_topology_selection.log",
    conda:
        "../envs/xsv.yaml"
    shell:
        "( xsv cat rows "
        '    --delimiter "\\t" '
        "    {input.rf_dists} | "
        "  xsv fmt --out-delimiter '\\t' "
        "  >{output.tsv} "
        ") 2>{log}"


rule join_raxml_ng_tsvs_per_individual_across_metrics:
    input:
        bsconverge="results/trees/{individual}.bsconverge.raxml_ng.max_missing_stable_topology_selection.tsv",
        rf_dist="results/trees/{individual}.rf_dist.max_missing_stable_topology_selection.tsv",
    output:
        tsv="results/trees/{individual}.max_missing_stable_topology_selection.tsv",
    log:
        "logs/trees/{individual}.max_missing_stable_topology_selection.join_metrics.log",
    conda:
        "../envs/xsv.yaml"
    shell:
        "( xsv join "
        '    --delimiter "\\t" '
        "    --left "
        "    individual,model,max_missing,tree_type,total_trees {input.rf_dist} "
        "    individual,model,max_missing,tree_type,total_trees {input.bsconverge} | "
        "  xsv select '!individual[1],model[1],max_missing[1],tree_type[1],total_trees[1]' | "
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
        support="results/gotree/{individual}/results/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.support.{tree_type}.raxml.collapsed.support",
        samples=config["samples"],
    output:
        support="results/trees/{individual}/{model}/max_{n_missing_cells}_missing/{individual}.{model}.max_{n_missing_cells}_missing.support.{tree_type}.collapsed.pdf",
    log:
        "logs/trees/{model}/max_{n_missing_cells}_missing/{individual}.support.{tree_type}.collapsed.log",
    conda:
        "../envs/ggtree.yaml"
    script:
        "../scripts/plot_support_tree.R"


rule plot_values_across_missingness:
    input:
        support_trees=expand(
            "results/gotree/{{individual}}/results/{model}/max_{n_missing_cells}_missing/{{individual}}.{model}.max_{n_missing_cells}_missing.support.{tree_type}.raxml.collapsed.support",
            model=lookup("raxml_ng/models", within=config),
            n_missing_cells=lookup(dpath="raxml_ng/max_missing", within=config),
            tree_type=["bestTree", "consensusTree"]
        ),
        cluster_info_dist=expand(
            "results/gotree/{{individual}}/results/{model}/max_{n_missing_cells}_missing/{{individual}}.{model}.max_{n_missing_cells}_missing.{type}.all.cluster_info_dist.tsv.gz",
            model=lookup("raxml_ng/models", within=config),
            n_missing_cells=lookup(dpath="raxml_ng/max_missing", within=config),
            type=["startTree", "mlTrees", "bootstraps"],
        )
    output:
        support_plot="results/trees/{individual}.tree_values_across_missingness.pdf",
    log:
        "logs/trees/{individual}.tree_values_across_missingness_plot.log",
    conda:
        "../envs/ggtree.yaml"
    resources:
        mem_mb=lambda wc, input: input.size_mb * 5
    script:
        "../scripts/plot_tree_values_across_missingness.R"


rule plot_support_values_and_branch_lengths:
    input:
        support_trees=expand(
            "results/gotree/{{individual}}/results/{model}/max_{n_missing_cells}_missing/{{individual}}.{model}.max_{n_missing_cells}_missing.support.{{tree_type}}.raxml.collapsed.support",
            model=lookup("raxml_ng/models", within=config),
            n_missing_cells=lookup(dpath="raxml_ng/max_missing", within=config),
        ),
        tsv=expand(
            "results/raxml_ng/{{individual}}/results/{model}/max_{n_missing_cells}_missing/{{individual}}.{model}.max_{n_missing_cells}_missing.bootstraps.raxml.bsconverge.tsv",
            model=lookup("raxml_ng/models", within=config),
            n_missing_cells=lookup(dpath="raxml_ng/max_missing", within=config),
        ),
    output:
        support_hist="results/trees/{individual}.{tree_type}.support.histogram.pdf",
        branch_length_hist="results/trees/{individual}.{tree_type}.branch_length.histogram.pdf",
        branch_length_ecdf="results/trees/{individual}.{tree_type}.branch_length.ecdf.pdf",
        data_plot="results/trees/{individual}.{tree_type}.support_vs_branch_length.full_data.pdf",
        summary_plot="results/trees/{individual}.{tree_type}.support_vs_branch_length.summary.pdf",
    log:
        "logs/trees/{individual}.{tree_type}.support_values_vs_branch_lengths.log",
    conda:
        "../envs/ggtree.yaml"
    script:
        "../scripts/plot_bootstrap_support_values_vs_branch_lengths.R"


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
