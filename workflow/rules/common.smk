import pandas as pd
import csv
import re
from os import path
from tempfile import TemporaryDirectory

from snakemake.utils import validate

#### validate config (read in in Snakefile)

validate(config, schema="../schemas/config.schema.yaml")


#### read in and validate samples.tsv

samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype={
            "sample_name": str,
            "individual": str,
            "sample_type": str,
            "platform": str,
        },
        comment="#",
    )
    .set_index("sample_name", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")


#### read in and validate units.tsv

units = (
    pd.read_csv(
        config["units"],
        sep="\t",
        dtype={
            "sample_name": str,
            "unit_name": str,
            "fq1": str,
            "cutadapt": str,
        },
        comment="#",
    )
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)

validate(units, schema="../schemas/units.schema.yaml")


#### constrain wildcards to enable unambiguous matching


wildcard_constraints:
    individual="|".join(samples["individual"].unique()),
    sample="|".join(samples["sample_name"]),
    sc="|".join(samples["sample_name"]),
    genotype="|".join(["hom_ref", "het", "hom_alt"]),
    ref_alt="[ACGT]_[ACGT]",


#### compile wanted workflow outputs


def get_final_output():
    final_output = []

    for individual in pd.unique(samples.individual):
        final_output.extend(
            expand(
                "results/final-calls/{ind}/{sc}.merged_bulk.prosolo.sorted.bcf",
                ind=individual,
                sc=samples.loc[
                    (samples["individual"] == individual)
                    & (samples["sample_type"] == "single_cell")
                ]["sample_name"],
            )
        )
    final_output.append("results/qc/multiqc.html")


#
#    for sample, unit in units.index:
#        row = units.loc[sample].loc[unit]
#        final_output.extend(
#            expand(
#                "results/qc/falco/{file}",
#                file=[
#                    row.fq1.replace(".fastq.gz", "/fastqc_report.html"),
#                    row.fq2.replace(".fastq.gz", "/fastqc_report.html"),
#                ],
#            )
#        )
#        final_output.extend(
#            expand(
#                "results/qc/falco/results/trimmed/{s}.{u}.{r}/fastqc_report.html",
#                s=sample,
#                u=unit,
#                r=["1", "2"],
#            )
#        )

    return final_output

#### helper functions

def get_single_cells_for_individual(individual):
    return samples.loc[ (samples["individual"] == individual) & (samples["sample_type"] == "single_cell"), "sample_name"]


#### input functions


def get_multiqc_input(wildcards):
    multiqc_input = []
    multiqc_input.extend(
        expand("results/qc/markdup/{sample}.metrics.txt", sample=samples.sample_name)
    )

    return multiqc_input


def get_sample_unit_fastqs(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    return [unit.fq1, unit.fq2]


def get_univec_reference_input(wildcards):
    if config["univec"]["activate"]:
        core = "_Core" if config["univec"]["core"] else ""
        return f"resources/reference/UniVec_Core/UniVec{core}.fa"
    else:
        return "/dev/null"


def get_merge_fastqs_input(wc):
    return expand(
        "results/trimmed/{sample}.{unit}.{read}.fastq.gz",
        sample=wc.sample,
        unit=units.loc[wc.sample, "unit_name"],
        read=wc.read,
    )


def get_processed_consensus_input(wildcards):
    if wildcards.read_type == "se":
        return "results/consensus/fastq/{}.se.fq".format(wildcards.sample)
    return [
        "results/consensus/fastq/{}.1.fq".format(wildcards.sample),
        "results/consensus/fastq/{}.2.fq".format(wildcards.sample),
    ]


def get_individual_samples(individual):
    return samples.loc[samples["individual"] == individual]["sample_name"]

def get_individual_samples_sample_type(individual, sample_type):
    return samples.loc[
        (samples["individual"] == individual)
        & (samples["sample_type"] == sample_type)
    ]["sample_name"]

def get_individual_bulk_samples_bam(wildcards):
    return expand("results/recal/{b}.sorted.bam", b=get_individual_samples_sample_type(wildcards.individual, "bulk"))


def aggregate_freebayes_region_calls_input(ext=".bcf"):
    def inner(wildcards):
        # decision based on content of output file
        # Important: use the method open() of the returned file!
        # This way, Snakemake is able to automatically download the file if it is generated in
        # a cloud environment without a shared filesystem.
        with checkpoints.create_freebayes_regions.get(
            individual=wildcards.individual
        ).output[0].open() as f:
            return expand(
                "results/candidate_calls/{ind}/{chr_and_region}.freebayes.norm{ext}",
                ind=wildcards.individual,
                chr_and_region=[chr_and_region.strip() for chr_and_region in f],
                ext=ext,
            )

    return inner


def get_all_scelestial_gts_for_individual(wildcards):
    single_cells = samples.loc[
        (samples["individual"] == wildcards.individual)
        & (samples["sample_type"] == "single_cell"),
        "sample_name",
    ]
    return expand(
        "results/scelestial/{individual}/{sc}.all_genotypes.tsv",
        individual=wildcards.individual,
        sc=single_cells,
    )


#### params functions


def get_cutadapt_parameters(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    try:
        cutadapt = unit["cutadapt"]
        if isinstance(cutadapt, str):
            return cutadapt
        return ""
    except KeyError:
        return ""


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample, platform=samples.loc[wildcards.sample, "platform"]
    )


def get_raxml_ng_mem_mb(wildcards, input):
    with open(input.log) as log:
        for line in log:
            if line.startswith("* Estimated memory requirements                : "):
                m = re.search(
                    r"\* Estimated memory requirements                : (\d+) MB", line
                )
                return int(m.group(1)) * 32


def get_raxml_ng_prefix(wildcards, output):
    return path.join(
        path.dirname(output[0]), path.basename(output[0]).rsplit(".raxml.")[0]
    )


def get_raxml_ng_invariant_sites(wildcards, input):
    with open(input.invariant_sites) as invariant_sites:
        num_sites=invariant_sites.readline().strip()
        if invariant_sites.readline():
            raise AssertionError(f"File '{input.invariant_sites}' should contain one value in one row. More than one row found.")
        return f"+ASC_FELS{{{num_sites}}}"


def get_raxml_ng_threads(wildcards, input):
    with open(input.log) as log:
        for line in log:
            if line.startswith("* Recommended number of threads / MPI processes: "):
                m = re.search(
                    r"\* Recommended number of threads / MPI processes: (\d+)", line
                )
                return int(m.group(1))
