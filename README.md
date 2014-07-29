ChIPseq-Analysis-Pipeline
=========================

This pipeline is written in Snakemake and takes SRA files and outputs read counts per input region.

####This pipeline is more of an example of how Snakemake can be used to quickly build pipelines than an actual off-the shelf pipeline.
The usage of this particular pipeline is quite specific to some of the work I'm currently doing.
The submit.sh is an example of how snakemake can be run on clusters (in my case, the NIH Biowulf clusters).

**One important aspect to remember with Snakemake is that the first rule ("all") is the last rule to be implemented.  This means that this first rule should include the output that you're expecting and Snakemake will attempt to create ONLY those outputs.  This is one way of controling which parts of the pipeline you want to implement and is a great for testing.**

Nonetheless, here's a quick description:

The pipeline can perform the following tasks:

	1) Convert SRA file to fastq file
	2) Generate QC output using fastqc
	3) Align reads using bowtie
	4) Run flagstat to calculate the number of reads
	5) Index the bam file
	6) Retrieve read counts per region (using htseq-count)

**Note**: You'll need the following software installed:

* [samtools](http://samtools.sourceforge.net/)
* [sratoolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)
* [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
* [bamtools](https://github.com/pezmaster31/bamtools)
* [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)

To run, all you need to do is

```
	qsub -l nodes=1:gpfs submit.sh
```
