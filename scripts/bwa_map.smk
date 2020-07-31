rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
        #"data/samples/{sample}.fastq.gz"
    output:
        temp("mapped_reads/{sample}.bam")
    threads: 2
    params:
        rg = r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    shell:
        "bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -h - > {output} 2> {log}" #注意{params.rg}两边的引号

