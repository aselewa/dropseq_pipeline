##Preparation of anotation files for the Drop-seq Pipelines



### 1. Obtain reference genomes



#### hg38
```bash
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz .
#combine into one single drop alternative haplotypes (_alt)
tar xvzf hg38.chromFa.tar.gz
cd chroms
rm *_alt.fa
cat *.fa > ../hg38_UCSC/hg38_ucsc.fa
rm -r chroms
```



#### mm10

```bash
#chromFa.tar.gz - The assembly sequence in one file per chromosome.
#    Repeats from RepeatMasker and Tandem Repeats Finder (with period
#    of 12 or less) are shown in lower case; non-repeating sequence is
#    shown in upper case.

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz .
tar xvzf chromFa.tar.gz
cat *.fa > mm10_UCSC/mm10_ucsc.fa
rm *.fa

```

#### danRer10 (zebrafish)

```bash
#Danio ri.. zebrafish

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/danRer10/bigZips/danRer10.fa.gz .

gunzip danRer10.fa.gz
```
The model system expresses additional markers that have to be included in the reference genome:

EGFP
>EGFP
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGGGCGGTGGAAGATCTGGGAATTCAAGGCCTCTCGAGCCTCTAGATTCTGCAGCCCTATAGC

mCherry
>mCherry
ATGGTGAGCAAGGGCGAGGAGGACAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCCGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGCTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAGCGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG


mouse c-myc
>m_cmyc
ATGCCCCTCAACGTGAACTTCACCAACAGGAACTATGACCTCGACTACGACTCCGTACAGCCCTATTTCATCTGCGACGAGGAAGAGAATTTCTATCACCAGCAACAGCAGAGCGAGCTGCAGCCGCCCGCGCCCAGTGAGGATATCTGGAAGAAATTCGAGCTGCTTCCCACCCCGCCCCTGTCCCCGAGCCGCCGCTCCGGGCTCTGCTCTCCATCCTATGTTGCGGTCGCTACGTCCTTCTCCCCAAGGGAAGACGATGACGGCGGCGGTGGCAACTTCTCCACCGCCGATCAGCTGGAGATGATGACCGAGTTACTTGGAGGAGACATGGTGAACCAGAGCTTCATCTGCGATCCTGACGACGAGACCTTCATCAAGAACATCATCATCCAGGACTGTATGTGGAGCGGTTTCTCAGCCGCTGCCAAGCTGGTCTCGGAGAAGCTGGCCTCCTACCAGGCTGCGCGCAAAGACAGCACCAGCCTGAGCCCCGCCCGCGGGCACAGCGTCTGCTCCACCTCCAGCCTGTACCTGCAGGACCTCACCGCCGCCGCGTCCGAGTGCATTGACCCCTCAGTGGTCTTTCCCTACCCGCTCAACGACAGCAGCTCGCCCAAATCCTGTACCTCGTCCGATTCCACGGCCTTCTCTCCTTCCTCGGACTCGCTGCTGTCCTCCGAGTCCTCCCCACGGGCCAGCCCTGAGCCCCTAGTGCTGCATGAGGAGACACCGCCCACCACCAGCAGCGACTCTGAAGAAGAGCAAGAAGATGAGGAAGAAATTGATGTGGTGTCTGTGGAGAAGAGGCAAACCCCTGCCAAGAGGTCGGAGTCGGGCTCATCTCCATCCCGAGGCCACAGCAAACCTCCGCACAGCCCACTGGTCCTCAAGAGGTGCCACGTCTCCACTCACCAGCACAACTACGCCGCACCCCCCTCCACAAGGAAGGACTATCCAGCTGCCAAGAGGGCCAAGTTGGACAGTGGCAGGGTCCTGAAGCAGATCAGCAACAACCGCAAGTGCTCCAGCCCCAGGTCCTCAGACACGGAGGAAAACGACAAGAGGCGGACACACAACGTCTTGGAACGTCAGAGGAGGAACGAGCTGAAGCGCAGCTTTTTTGCCCTGCGTGACCAGATCCCTGAATTGGAAAACAACGAAAAGGCCCCCAAGGTAGTGATCCTCAAAAAAGCCACCGCCTACATCCTGTCCATTCAAGCAGACGAGCACAAGCTCACCTCTGAAAAGGACTTATTGAGGAAACGACGAGAACAGTTGAAACACAAACTCGAACAGCTTCGAAACTCTGGTGCATAA



```bash
cat danRer10.fa EGFP.fa mCherry.fa c_myc.fa > danRer10_RG_cmyc.fa

```


#### panTro5

```bash


mkdir panTro5_UCSC
cd panTro5_UCSC
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/panTro5/bigZips/panTro5.fa.gz .
tar -xvzf panTro5.fa.gz

panTro5.fa.gz  


```



### 2. Obtain transcript annotation

hg38
```bash
Release 27 (GRCh38.p10)
ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
#also prepare a gtf file version without pseudogenes

less gencode.v27.annotation.gtf |grep -v "pseudogene" > gencode.v27.annotation.no.pseudo.gtf


```


mm10
```bash
Release M15 (GRCm38.p5)
ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M15/gencode.vM15.annotation.gtf.gz



```

danRer10
```bash
# downloaded danRer10 from UCSC
RefGene_danRer10-2017-08-28.gtf

#in addition add the annotations for the reporter genes EGFP and mCherry, and the overexpressed oncogene c-myc
#lengths
EGFP	danRer10_refGene	CDS	1	778	.	+	.	gene_id	"EGFP";	transcript_id "EGFP";
EGFP	danRer10_refGene	exon	1	778	.	+	.	gene_id "EGFP"; transcript_id "EGFP";
mCherry	danRer10_refGene	CDS	1	708	.	+	.	gene_id "mCherry"; transcript_id "mCherry";
mCherry danRer10_refGene	exon	1	708	.	+	.	gene_id "mCherry"; transcript_id "mCherry";
m_cmyc  danRer10_refGene  CDS 1 1320  . + . gene_id "m_cmyc"; transcript_id "m_cmyc";
m_cmyc  danRer10_refGene  exon 1 1320  . + . gene_id "m_cmyc"; transcript_id "m_cmyc";



#just add the below:
EGFP    AddedGene        exon    1       778     .       +       0       gene_id "EGFP"; transcript_id "EGFP";
mCherry AddedGene        exon    1       708     .       +       0       gene_id "mCherry"; transcript_id "mCherry";
cMyc  AddedGene        exon    1       1320    .       +       0       gene_id "m_cmyc"; transcript_id "m_cmyc";


#combine into additional_sequences.gtf

cat RefGene_danRer10-2017-08-28.gtf	 additional_sequences.gtf > RefGene_danRer10-2017-08-28_RG_cmyc.gtf



```



### 3. Prepare STAR indeces

Read length is 62 bp for the current DropSeq set up. If lengths of reads changes (reads in R2 file), this needs to be changed with the appropriate length


#### hg38
```bash
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir hg38_noalt_juncGencodeV27_61 --genomeFastaFiles hg38_UCSC/hg38_ucsc.fa --sjdbGTFfile hg38_UCSC/gencode.v27.annotation.gtf --sjdbOverhang 61


```

#### mm10

```bash
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir mm10_noalt_juncGencodemV15_61 --genomeFastaFiles mm10_UCSC/mm10_ucsc.fa --sjdbGTFfile mm10_UCSC/gencode.vM15.annotation.gtf --sjdbOverhang 61


```
#### Zebrafish (danRer10_RGcMyc)
This index includes the fluorochromes EGFP and mCherry as well as the oncogene c-myc (from mouse).

STAR --runThreadN 2 --runMode genomeGenerate --genomeDir dr10_noalt_juncRefGene_RG_myc_61 --genomeFastaFiles danRer10_UCSC_RG/danRer10_RG_cmyc.fa --sjdbGTFfile danRer10_UCSC_RG/RefGene_danRer10-2017-08-28_RG_cmyc.gtf --sjdbOverhang 61
```

### 4. Obtain files for specific features
Get the interval files for the following features:
-- genes
-- exons
-- introns
-- non-genic
-- ribosomal


#### hg38




just use bedtools coverage for this
NO: use CollectRnaSeqMetrics


for riboomal list: Represents a list of intervals against a reference sequence that can be written to and read from a file. The file format is relatively simple and reflects the SAM alignment format to a degree. A SAM style header must be present in the file which lists the sequence records against which the intervals are described. After the header the file then contains records one per line in text format with the following values tab-separated: Sequence name, Start position (1-based), End position (1-based, end inclusive), Strand (either + or -), Interval name (an, ideally unique, name for the interval),


```bash
#prepare a flat ref file
 gtfToGenePred -genePredEx gencode.v27.annotation.gtf gencode.v27.annotation.refFlat
 less gencode.v27.annotation.refFlat |awk '{print($12"\t"$0)}' | cut -f1-11 >gencode.v27.annotation.NameCol.refFlat

 less gencode.v27.annotation.gtf |grep  "gene_type \"rRNA\"" |head



samtools faidx hg38_ucsc.fa
samtools view -ht hg38_ucsc.fa.fai hg38_rRNA.txt > header #this fails to make a proper file, but yields the header
cat header hg38_rRNA.txt > hg38_ribosome.interval_list





Gene_regions.bed			gencode.v27.annotation.gtf		header					hg38_ribosome.interval_list		hg38_ucsc.fa.fai

java -jar picard.jar CollectRnaSeqMetrics \
      I=input.bam \
      O=output.RNA_Metrics \
      REF_FLAT=ref_flat.txt \
      STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
      RIBOSOMAL_INTERVALS=ribosomal.interval_list


      java -jar picard CollectRnaSeqMetrics \
            I=/Volumes/Data2/Drop_seq_pipeline/assigned_sorted.bam \
            O=output.RNA_Metrics \
            REF_FLAT=STAR_indeces/hg38_UCSC/gencode.v27.annotation.NameCol.refFlat \
            STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
            RIBOSOMAL_INTERVALS=STAR_indeces/hg38_UCSC/hg38_ribosome.interval_list

```

notes:  "Picard version 1.97 uses a wrapper script symlinked to bin/: picard.  The wrapper script takes as arguments the particular jar file you want to run followed by any arguments for that jar.  For example, 'picard ViewSam.jar --help'."

java -Xmx12G picard CollectRnaSeqMetrics \
      I=/Volumes/Data2/Drop_seq_pipeline/assigned_sorted.bam \
      O=output.RNA_Metrics \
      REF_FLAT=STAR_indeces/hg38_UCSC/gencode.v27.annotation.NameCol.refFlat \
      STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
      RIBOSOMAL_INTERVALS=STAR_indeces/hg38_UCSC/hg38_ribosome.interval_list



 --help
