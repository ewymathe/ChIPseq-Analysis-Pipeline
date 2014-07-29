import config
import glob
import os
import gzip

def get_read_len(wildcards):
     f=gzip.open((wildcards + "_1.fastq.gz"),'rb')
     seq=f.readline()
     seq=f.readline()
     f.close()
     return(str(len(seq)-1))

def getpath(wildcards):
     mysplit=wildcards.split('/')
     return(mysplit[0]+"/"+mysplit[1])

# The input samples need to be described here:
SAMPLES = [fname.split('.')[0] for fname in glob.glob('*/*/*.bam')]
print(SAMPLES)

# This rule, the first rule, is the last rule run in Snakemake.  
# So it should include the outputs that you're expecting and it will only create these outputs.
rule all:
     params: batch='-l nodes=1:8'
     input: expand('{sample}.counts', sample=SAMPLES)

rule fastqdump:
    """Run fastq-dump to convert sra format to fastq files, creating .fq"""
    input: "{sample}.sra"
    output: "{sample}_1.fastq.gz"
    params: batch = "-l nodes=1:g8", mydir = lambda wildcards: getpath(wildcards.sample)
    threads: 1
    shell: """module load sratoolkit; fastq-dump -split-files --gzip {input} -O {params.mydir} """

rule runbowtie:
    """Run alignment using bowtie, creating .bam"""
    input: "{sample}_1.fastq.gz"
    output: "{sample}.bam"
    params: batch = "-l nodes=1:gpfs",
            readlen = lambda wildcards: get_read_len(wildcards.sample),
            bowtie_params = "--sam --best --all --strata -m1 -n2 --threads=16",
            bowtie_reference = config.BOWTIE_REFERENCE
    threads: 16
    shell: """module load bowtie; gunzip -c {input} | bowtie {params.bowtie_params} -l{params.readlen} {params.bowtie_reference} - | \
       samtools view -Sb -F4 - 2> {wildcards.sample}.log > {wildcards.sample}unsorted.bam && \
       samtools sort {wildcards.sample}unsorted.bam {output} && mv {output}.bam {output} && \
       rm {wildcards.sample}unsorted.bam"""

rule getstats:
     """ Run samtools flagstat to the number of aligned reads """
     input: "{sample}.bam"
     output: "{sample}.bam.stat"
     params: batch = "-l nodes=1:g4"
     threads: 1
     shell: """samtools flagstat {input} > {output}"""

rule index:
     """ Run bamtools index to index bam files """
     input: "{sample}.bam"
     output: "{sample}.bam.bai"
     params: batch = "-l nodes=1:gpfs"
     threads: 1
     shell: "module load bamtools; bamtools index -in {input}"

rule getcounts:
     """ Run htseq-count to calculate the number counts that overlap with input regions """
     input: "{sample}.bam"
     output: "{sample}.counts" 
     params: batch = "-l nodes=1:gpfs"
     threads: 1
     shell: """ htseq-count {input} myregions.gtf > {output} """ 
