configfile: "config.yaml"

SAMPLES = config["samples"]
THREADS = 16

rule all:
    input:
        expand("qc_reports/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        expand("qc_reports/{sample}_R2_001_fastqc.html", sample=SAMPLES),
        expand("kallisto_out/{sample}/counts.mtx", sample=SAMPLES),
        "busco_results/short_summary_busco_transcriptome.txt"

rule kallisto_index:
    input:
        config["transcriptome"]
    output:
        "transcriptome.idx"
    conda:
        "envs/kallisto.yaml"
    threads: THREADS
    shell:
        "kallisto index -i {output} {input}"

rule transcripts_to_genes:
    input:
        config["transcriptome"]
    output:
        "transcripts_to_genes.txt"
    shell:
        "grep '^>' {input} | sed 's/>//' | awk '{{print $1\t$1}}' > {output}"

rule kallisto_bus:
    input:
        idx="transcriptome.idx",
        r1="reads/{sample}_R1_001.fastq.gz",
        r2="reads/{sample}_R2_001.fastq.gz"
    output:
        bus="kallisto_out/{sample}/output.bus"
    params:
        outdir="kallisto_out/{sample}",
        tech="10xv2"
    conda:
        "envs/kallisto.yaml"
    threads: THREADS
    shell:
        """
        mkdir -p {params.outdir} && \
        kallisto bus -i {input.idx} -o {params.outdir} -x {params.tech} -t {threads} {input.r1} {input.r2}
        """

rule bustools_pipeline:
    input:
        bus="kallisto_out/{sample}/output.bus",
        t2g="transcripts_to_genes.txt",
        whitelist=config["whitelist"]
    output:
        mtx="kallisto_out/{sample}/counts.mtx",
        genes="kallisto_out/{sample}/counts.genes.txt",
        barcodes="kallisto_out/{sample}/counts.barcodes.txt"
    conda:
        "envs/bustools.yaml"
    threads: THREADS
    shell:
        """
        cd kallisto_out/{wildcards.sample} && \
        bustools correct -w ../../{input.whitelist} -o output.corrected.bus output.bus && \
        bustools sort -t {threads} -o output.sorted.bus output.corrected.bus && \
        bustools count -o counts -g ../../{input.t2g} -e matrix.ec -t transcripts.txt --genecounts output.sorted.bus
        """

rule fastqc:
    input:
        r1="reads/{sample}_R1_001.fastq.gz",
        r2="reads/{sample}_R2_001.fastq.gz"
    output:
        "qc_reports/{sample}_R1_001.fastqc.html",
        "qc_reports/{sample}_R2_001.fastqc.html"
    conda:
        "envs/fastqc.yaml"
    threads: THREADS
    shell:
        "fastqc {input.r1} {input.r2} --outdir=qc_reports"

rule run_busco:
    input:
        "transcriptome.fa"
    output:
        "busco_results/short_summary_busco_transcriptome.txt"
    conda:
        "envs/busco.yaml"
    threads: THREADS
    shell:
        """
        mkdir -p busco_results && \
        busco -i {input} -o busco_transcriptome -m transcriptome -l metazoa_odb10 -f --cpu 16 && \
        cp $(ls busco_transcriptome/short_summary*.txt | head -n 1) {output}
        """
