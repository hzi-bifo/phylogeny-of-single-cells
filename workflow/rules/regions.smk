rule build_sample_regions:
    input:
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bai",
    output:
        "results/regions/{individual}/{sample}.mosdepth.global.dist.txt",
        "results/regions/{individual}/{sample}.quantized.bed.gz",
        summary="results/regions/{individual}/{sample}.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth/regions/{individual}_{sample}.log",
    params:
        extra="--no-per-base",
        quantize="1:",
    resources:
        runtime=lambda wc, attempt: attempt * 59 - 1,
    wrapper:
        "v1.21.1/bio/mosdepth"


rule merge_individual_regions:
    input:
        lambda wc: expand(
            "results/regions/{{individual}}/{sample}.quantized.bed.gz",
            sample=get_individual_samples(wc.individual),
        ),
    output:
        "results/regions/{individual}.covered_regions.bed",
    log:
        "logs/regions/{individual}_covered_regions.log",
    conda:
        "../envs/bedtools.yaml"
    resources:
        runtime=lambda wc, attempt: attempt * 59 * 2 - 1,
    shell:
        "zcat {input} | sort -k1,1 -k2,2n - | mergeBed -i - -d 15000 > {output} 2> {log}"


rule filter_individual_regions:
    input:
        covered="results/regions/{individual}.covered_regions.bed",
        fai=rules.full_reference_faidx.output,
    output:
        "results/regions/{individual}.covered_regions.filtered.bed",
    params:
        chroms=config["ref"].get("n_chromosomes", 25),
    log:
        "logs/regions/{individual}.covered_regions.filtered.log",
    shell:
        "cat {input.covered} | grep -f <(head -n {params.chroms} {input.fai} | "
        'awk \'{{print "^"$1"\\t"}}\') '
        "> {output} 2> {log}"

#### alternative minimum coverage regions determination introduced during
#### fine-tuning of raxml-ng input
# TODO: use the rule path below instead of build_sample_regions above for the
#       the next workflow release after the publication

rule build_sample_regions_cov:
    input:
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bai",
    output:
        "results/regions/{individual}/{sample}.min_cov_for_candidate_sites.mosdepth.global.dist.txt",
        "results/regions/{individual}/{sample}.min_cov_for_candidate_sites.quantized.bed.gz",
        summary="results/regions/{individual}/{sample}.min_cov_for_candidate_sites.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth/regions/{individual}_{sample}.min_cov_for_candidate_sites.log",
    params:
        extra="--no-per-base",
        quantize=f"{config['min_cov_for_candidate_sites']}:",
    resources:
        runtime=lambda wc, attempt: attempt * 59 - 1,
    wrapper:
        "v1.21.1/bio/mosdepth"


rule filter_individual_regions_cov:
    input:
        covered="results/regions/{individual}/{sample}.min_cov_for_candidate_sites.quantized.bed.gz",
        fai=rules.full_reference_faidx.output,
    output:
        "results/regions/{individual}/{sample}.min_cov_for_candidate_sites.covered_regions.filtered.bed",
    params:
        chroms=config["ref"].get("n_chromosomes", 25),
    log:
        "logs/regions/{individual}.{sample}.covered_regions.filtered.log",
    shell:
        "zcat {input.covered} | grep -f "
        "  <(head -n {params.chroms} {input.fai} | awk '{{print \"^\"$1\"\\t\"}}') "
        "> {output} 2> {log}"


rule merge_individual_cell_regions_cov:
    input:
        beds=lambda wc: expand(
            "results/regions/{{individual}}/{sample}.min_cov_for_candidate_sites.covered_regions.filtered.bed",
            sample=get_individual_samples_sample_type(wc.individual, "single_cell"),
        ),
        fai=rules.full_reference_faidx.output,
    output:
        "results/regions/{individual}.min_cov_for_candidate_sites.covered_regions.single_cells.bed",
    log:
        "logs/regions/{individual}.min_cov_for_candidate_sites.covered_regions.sindle_cells.log",
    params:
        min_cells=lambda wc, input: max(round(len(input.beds)/2), 5),
    conda:
        "../envs/bedtools.yaml"
    resources:
        runtime=lambda wc, attempt: 59 * 2 * attempt,
    shell:
        "(cat {input.beds} | "
        " sort -k1,1 -k2,2n - | "
        " bedtools genomecov -i - -g {input.fai} -bg -max {params.min_cells} | "
        " awk '$4 >= {params.min_cells}' "
        " >{output} "
        ") 2>{log} "


rule merge_individual_bulk_regions_cov:
    input:
        beds=lambda wc: expand(
            "results/regions/{{individual}}/{sample}.min_cov_for_candidate_sites.covered_regions.filtered.bed",
            sample=get_individual_samples_sample_type(wc.individual, "bulk"),
        ),
        fai=rules.full_reference_faidx.output,
    output:
        "results/regions/{individual}.min_cov_for_candidate_sites.covered_regions.bulks.bed",
    log:
        "logs/regions/{individual}.min_cov_for_candidate_sites.covered_regions.bulks.log",
    params:
        n_samples=lambda wc, input: len(input.beds),
    conda:
        "../envs/bedtools.yaml"
    resources:
        runtime=lambda wc, attempt: 59 * 2 * attempt,
    shell:
        "(cat {input.beds} | "
        " sort -k1,1 -k2,2n - | "
        " bedtools genomecov -i - -g {input.fai} -bg -max {params.n_samples} | "
        " awk '$4 >= {params.n_samples}' "
        " >{output} "
        ") 2>{log} "


rule consolidate_candidate_regions:
    input:
        "results/regions/{individual}.min_cov_for_candidate_sites.covered_regions.bed",
    output:
        "results/regions/{individual}.min_cov_for_candidate_sites.covered_regions.consolidated.bed",
    log:
        "logs/regions/{individual}.min_cov_for_candidate_sites.covered_regions.consolidated.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools merge -i {input} -d 15000 "
        ">{output} 2>{log}"


rule get_total_individual_cov:
    input:
        single_cells="results/regions/{individual}.min_cov_for_candidate_sites.covered_regions.single_cells.bed",
        bulks="results/regions/{individual}.min_cov_for_candidate_sites.covered_regions.bulks.bed",
        fai=rules.full_reference_faidx.output,
    output:
        "results/regions/{individual}.min_cov_for_candidate_sites.covered_sites.txt",
    log:
        "logs/regions/{individual}.min_cov_for_candidate_sites.covered_sites.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "(cat {input.single_cells} {input.bulks} | "
        " sort -k1,1 -k2,2n - | "
        " bedtools genomecov -i - -g {input.fai} -bg -max 2 | "
        " awk '$4 >= 2' | "
        " awk '{{sum += ($3 - $2)}} END {{print sum}}' "
        " >{output} "
        ") 2>{log} "
