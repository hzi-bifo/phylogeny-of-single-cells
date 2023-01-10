import tempfile
from snakemake.shell import shell

# Extract arguments
extra = snakemake.params.get("extra", "")
reference = snakemake.input.get("ref")
embed_ref = snakemake.params.get("embed_ref", False)
java_opts = snakemake.params.get("java_opts", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=False)

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "(gatk --java-options '{java_opts}' ApplyBQSR"
        " --input {snakemake.input.bam}"
        " --bqsr-recal-file {snakemake.input.recal_table}"
        " --reference {reference}"
        " {extra}"
        " --tmp-dir {tmpdir}"
        " --output {output} ) {log}"
    )
