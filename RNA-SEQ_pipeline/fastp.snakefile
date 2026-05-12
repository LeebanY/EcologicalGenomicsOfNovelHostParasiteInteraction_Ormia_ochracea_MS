SAMPLES, = glob_wildcards("{sample}_1.fq.gz")
configfile: "kallisto.yml"
rule all:
        input:
                expand("cleaned/{sample}_1.cleaned.fastq", sample=SAMPLES)

ruleorder: fastp_pe > single_end_fastp

rule single_end_fastp:
        input:
                in_read1="{sample}_1.fq.gz"
        output:
                out_read1="cleaned/{sample}_1.cleaned.fastq"
        threads: 16
        log:
                "single_end_fastp.{sample}.log"
        shell:
                "fastp -i {input.in_read1} -o {output.out_read1}"

rule fastp_pe:
        input:
                in_read1="{sample}_1.fq.gz", in_read2="{sample}_2.fq.gz"
        output:
                out_read1="cleaned/{sample}_1.cleaned.fastq", out_read2="cleaned/{sample}_2.cleaned.fastq"
        threads: 16
        log:
                "fast_pe.{sample}.log"
        shell:
                "fastp -i {input.in_read1} -I {input.in_read2} -o {output.out_read1} -O {output.out_read2}"
