import pandas as pd
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


#### compile wanted workflow outputs


def get_final_output():
    final_output = []

    final_output.extend(
        expand(
            "results/candidate-calls/{individual}.freebayes.bcf",
            individual=pd.unique(samples.individual)
        )
    )
    final_output.append("results/qc/multiqc.html")
    for sample, unit in units.index:
        row = units.loc[sample].loc[unit]
        final_output.extend(
            expand(
                "results/qc/falco/{file}",
                file=[
                    row.fq1.replace(".fastq.gz", "/fastqc_report.html"),
                    row.fq2.replace(".fastq.gz", "/fastqc_report.html"),
                ],
            )
        )
        final_output.extend(
            expand(
                "results/qc/falco/results/trimmed/{s}.{u}.{r}/fastqc_report.html",
                s=sample,
                u=unit,
                r=["1", "2"],
            )
        )

    return final_output


#### input functions


def get_multiqc_input(wildcards):
    multiqc_input = []
    multiqc_input.extend(
        expand(
            "results/qc/markdup/{sample}.metrics.txt",
            sample=samples.sample_name
        )
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
