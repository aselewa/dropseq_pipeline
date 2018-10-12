## dropseq_pipeline

This Drop-seq pipeline is designed to process data from fastq files to a digital expression matrix (dge).
The pipeline is based on Snakemake and is currently designed to run on the uchicago rcc cluster **midway2**. Below is a description of how to set up the project folders and to start the analysis.
Each experiment differs and the pipeline might need to be adjusted to accommodate such individual differences.

Below are steps that are common to all experiments and outlines of analysis that are commonly used.

## Preparing Environment

#### 1. Create a compute environment using conda (This step can be skipped when re-running an analysis)

The environment needs to be created only once. It will be activated when running the dropseq pipeline.

```bash
module load Anaconda3

conda env create --file .../dropseq_pipeline/environment.yaml
```

To update the environment, you can run the following command:
```bash
conda env update --file .../dropseq_pipeline/environment.yaml
```

#### 2. Edit the config_hg38.yaml file

Edit the config_hg38.yaml file so that Snakemake can find the right files

## Prepare data:

#### 1. Create a project directory in your directory on midway2
```bash
mkdir your_project
```
#### 2. Create directory fastq in 'your_project' directory and add fastq files
```bash
cd your_project

mkdir data
cd data/

mkdir fastq
cd fastq/
#only include the fastq files included in a single run, both read 1 and read2
cp path/to/fastq/*fastq.gz .

cd ../../
```


## Run dropseq pipeline

This command will run the Submit_snakemake.sh and pass the location of your project directory.

```bash
.../dropseq_pipeline/snakemake.batch "--config proj_dir=/project2/PI/CNETID/Path/to/your/dir/"
```



