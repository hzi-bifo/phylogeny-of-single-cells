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


#### compile wanted workflow outputs


def get_final_output():
    final_output = []

    for sample, unit in units.index:
        final_output.extend(
            expand(
                "results/trimmed/{s}.{u}.{r}.fastq.gz", s=sample, u=unit, r=["1", "2"]
            )
        )

    final_output.append("results/qc/multiqc.html")
    final_output.append("resources/reference/full_reference.pac")

    return final_output


#### input functions


def get_sample_unit_fastqs(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    return [unit.fq1, unit.fq2]


def get_univec_reference_input(wildcards):
    if config["univec"]["activate"]:
        core = "_Core" if config["univec"]["core"] else ""
        return f"resources/reference/UniVec_Core/UniVec{core}.fa"
    else:
        return "/dev/null"


def get_multiqc_input(wildcards):
    multiqc_input = []
    for sample, unit in units.index:
        row = units.loc[sample].loc[unit]
        multiqc_input.extend(
            expand(
                "results/qc/fastqc/{file}",
                file=[
                    row.fq1.replace(".fastq.gz", "_fastqc.zip"),
                    row.fq2.replace(".fastq.gz", "_fastqc.zip"),
                ],
            )
        )
        multiqc_input.extend(
            expand(
                "results/qc/fastqc/results/trimmed/{s}.{u}.{r}_fastqc.zip",
                s=sample,
                u=unit,
                r=["1", "2"],
            )
        )

    return multiqc_input


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
