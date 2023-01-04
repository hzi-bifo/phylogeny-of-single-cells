import tempfile
from snakemake.shell import shell

log = snakemake.log_fmt_shell()

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "( picard MarkDuplicates"  # Tool and its subcommand
        "   {snakemake.params.java_opts}"  # Automatic java option
        "   {snakemake.params.extra}"  # User defined parmeters
        "   {snakemake.params.bams}"  # Input bam(s)
        "   --TMP_DIR {tmpdir}"
        "   --OUTPUT {snakemake.output.bam}"  # Output bam
        "   --METRICS_FILE {snakemake.output.metrics}"  # Output metrics
        ") {log} "  # Logging
    )