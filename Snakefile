#Snakefile for Dropseq analysis - this should be used first to try and infer number of cells in sample


import glob
import os

# Configuration ----------------------------------------------------------------

#expected number of cells (upper limit)
cell_num = config["cell_num"]

#cell barcode UMI structure
barcode = config["barcode"]

#genome_index
GenomeIndex = config["genome_index"]
#gene file
txn_file = config["txn_file"]

pd = config["proj_dir"]
data = "data/"
output =  "output/"
fastq_dir = data + "fastq/"
fastqc_dir = output + "fastqc/"
fastq_merged = data + "fastq_merged/"
fastq_extr = data + "fastq_extr/"
cell_stats = data + "cell_stats/"
aligned = data + "aligned/"
annotation_data = data + "annotation_data/"
sorted_reads = data + "sorted_reads/"
assigned = data + "assigned/"
dge_data = output + "dge_data/"
qc_data = output + "qc_data/"
scripts = pd + "Scripts/"
aggregate = data + "aggregate/"
code = "Code/"


##make sure the project directory actually exists
#sassert os.path.exists(pd), "Project directory exists"

# Directory to send log files. Needs to be created manually since it
# is not a file created by a Snakemake rule.
dir_log = config["dir_log"]
if not os.path.isdir(dir_log):
    os.mkdir(dir_log)
# input data # might need to be changed to be universal

samples = set(glob_wildcards(fastq_dir + "{samples}_R1_001.fastq.gz").samples)
percent = config["downsample"]
type = config["type"]

rule all:
    input:
        #cell_stats + "whitelist.txt",
        #fastq_extr + "combined_r2_extracted.fastq.gz",
        #aligned + "Aligned.SortedByCoordinate.out.bam",
        expand("data/aligned/Aligned.SortedByCoordinate.out.{type}.bam",type=type),
        expand("data/assigned/Aligned.SortedByCoordinate.out.{type}.bam.featureCounts.bam",type=type),
        #expand("data/assigned/Aligned.featureCounts.downsampled_{id}_{type}.bam",id=percent,type=type),
        expand("data/sorted_reads/assigned_sorted_{type}.bam",id=percent,type=type),
        expand("data/sorted_reads/assigned_sorted_{type}.bam.bai",id=percent,type=type),
        expand("output/dge_data/counts_{type}.tsv.gz",id=percent,type=type),
        #cell_stats + "whitelist.txt",
        #fastqc_dir + "combined_r1_fastqc.zip",
        #qc_data + "Nucleotide_frequency_in_UMIs.pdf",
        #qc_data + "Nucleotide_frequency_in_cell_barcode.pdf"


#fastqc will be run on merged files
rule fastqc:
    input:
        fastq_merged + "combined_r{read_num}.fastq.gz"
    output:
        fastqc_dir + "combined_r{read_num}_fastqc.html",
        fastqc_dir + "combined_r{read_num}_fastqc.zip"
    params:
        outdir = fastqc_dir
    shell:
        "fastqc -o {params.outdir} {input}"

rule unzip:
    input:
        fastqc_dir + "combined_r1_fastqc.zip"
    output:
        fastqc_dir + "combined_r1_fastqc/fastqc_data.txt"
    shell:
        "unzip {input}"

rule barcode_qc:
    input:
        fastqc_dir + "combined_r1_fastqc/fastqc_data.txt"
    output:
        qc_data + "Nucleotide_frequency_in_UMIs.pdf",
        qc_data + "Nucleotide_frequency_in_cell_barcode.pdf"
    shell:
        "shell Scripts/Calculate-nuc-freq.sh {input}"

rule fastq_merge_r1:
    input:
        expand(fastq_dir + "{sample}_R1_001.fastq.gz", sample = samples)
    output:
        fastq_merged + "combined_r1.fastq.gz"
    shell:
        "zcat {input} | gzip --stdout > {output}"

rule fastq_merge_r2:
    input:
        expand(fastq_dir + "{sample}_R2_001.fastq.gz", sample = samples)
    output:
        fastq_merged + "combined_r2.fastq.gz"
    shell:
        "zcat {input} | gzip --stdout  > {output}"


rule umi_create_whitelist:
    input:
        fastq_merged + "combined_r1.fastq.gz"
    output:
        cell_stats + "whitelist.txt"
    params:
        cell_num = cell_num,
        bc = barcode
    shell:
        "umi_tools whitelist --stdin {input} --bc-pattern={params.bc} --plot-prefix=Whitelist_stats --extract-method=string --expect-cells={params.cell_num} --log2stderr > {output}"

rule umi_extract_bc_and_umi:
    input:
        r1 = fastq_merged + "combined_r1.fastq.gz",
        r2 = fastq_merged + "combined_r2.fastq.gz",
        wl = cell_stats + "whitelist.txt"
    output:
        r1_ext = temp(fastq_extr + "combined_r1_extracted.fastq.gz"),
        r2_ext = fastq_extr + "combined_r2_extracted.fastq.gz"
    params:
        bc = barcode
    shell:
        "umi_tools extract --bc-pattern={params.bc} --stdin {input.r1} --stdout {output.r1_ext} --read2-in {input.r2} --read2-out={output.r2_ext}  --error-correct-cell --filter-cell-barcode --whitelist={input.wl}"


#allows for proportion of 0.1 mismatches (e.g. 10*0.1=1 mismatch in 10bp adaptor match, bp below zero are rounded to 0)
# sequences to match are poly at the 3 prime end (at least 6) and TSO oligo on the 5 prime end
rule trim_read2:
    input:
        r1 = fastq_extr + "combined_r2_extracted.fastq.gz"
    output:
        r1_trim = fastq_extr + "combined_r2_trimmed.fastq.gz",
    params:
        min_len = 30
    shell:
        "cutadapt --minimum-length {params.min_len} -a AAAAAA -g AAGCAGTGGTATCAACGCAGAGTGAATGGG -o {output} {input}"

rule align:
    input:
        fq = fastq_extr + "combined_r2_trimmed.fastq.gz",
        ref_genome = GenomeIndex
    output:
        aligned + "Aligned.SortedByCoordinate.out.bam"
    threads: 4
    shell:
        "STAR --runThreadN {threads} --genomeDir {input.ref_genome} --readFilesIn {input.fq} --readFilesCommand zcat --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > {output}"

rule duplicate_bams:
    input: 
       bam = aligned + "Aligned.SortedByCoordinate.out.bam"
    output:
       bam_out = aligned + "Aligned.SortedByCoordinate.out.{type}.bam"
    shell:
        "cp {input.bam} {output.bam_out}"

# assign reads to transcripts
rule reads_to_transcripts:
    input:
        bam = aligned + "Aligned.SortedByCoordinate.out.{type}.bam",
        features = txn_file
    output:
        assigned_feat = assigned + "gene_assigned.{type}.summary",
        bam_out = assigned + "Aligned.SortedByCoordinate.out.{type}.bam.featureCounts.bam"
    threads: 1
    shell:
        "featureCounts -a {input.features} -o {output.assigned_feat} -R BAM {input.bam} -T {threads} -t {wildcards.type} -g gene_id"

downsample reads
rule downsample:
    input:
        assigned + "Aligned.SortedByCoordinate.out.{type}.bam.featureCounts.bam"
    output:
        assigned  + "Aligned.featureCounts.downsampled_{id}_{type}.bam"
    shell:
        "samtools view -s {wildcards.id} -b {input} > {output}"

rule sort_bams:
    input:
        assigned + "Aligned.SortedByCoordinate.out.{type}.bam.featureCounts.bam"
    output:
        sorted_reads + "assigned_sorted_{type}.bam"
    shell:
        "samtools sort -o {output} -O bam {input} -T data/sorted_reads/{wildcards.type}_temp"

rule index_bams:
    input:
        sorted_reads + "assigned_sorted_{type}.bam"
    output:
        sorted_reads + "assigned_sorted_{type}.bam.bai"
    shell:
        "samtools index {input}"

rule make_DGE_matrix:
    input:
        sorted_reads + "assigned_sorted_{type}.bam",
        sorted_reads + "assigned_sorted_{type}.bam.bai"
    output:
        dge_data + "counts_{type}.tsv.gz"
    shell:
        "umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --per-cell -I {input} -S {output}"
