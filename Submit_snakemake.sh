#!/bin/bash

source activate dropseq

snakemake \
    -kp \
    --ri \
    -j 28 \
    --cluster-config /project2/xinhe/alan/dropseq_pipeline/cluster.json \
    -c "sbatch \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --tasks-per-node={cluster.tasks} \
        --partition=xinhe \
	--account=pi-xinhe \
        --job-name={cluster.name} \
	--output={cluster.logfile}" \
    $*
