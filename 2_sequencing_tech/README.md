# Sequencing Technology Comparison

## Requirements

* [Conda with Python 3 environment](https://conda.io/en/latest/miniconda.html)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/getting\_started/installation.html)

## Running

If running the program from within this directory, you can run the single command:

```
snakemake --use-conda
```

You may want to use the -j parameter to run on multiple cpu cores.

## Results

The outputs for the tables and figures are at:

* Table 4
 * results/filtering\_stats.txt
* Table 5
 * results/preprocessing\_stats.txt
* Table 6
 * results/completeness\_stats.txt 
* Table 7
 * results/shared\_complete\_transcripts/4gb\_percent.txt
* Figure 3
 * results/exome\_genome\_coverage.pdf
* Figure 4
 * results/gap\_vs\_covered\_4gb\_BGISEQ.pdf
 * results/gap\_vs\_covered\_4gb\_HiSeq.pdf
* Figure 5
 * results/gc\_bias/ERR1831364\_gc\_bias.pdf
 * results/gc\_bias/SRR1261168\_gc\_bias.pdf
 * results/gc\_bias/SRR950078\_gc\_bias.pdf
* Figure 6
 * results/kallisto\_plots/SRR1261168\_vs\_ERR1831364\_kallisto\_4gb.pdf
 * results/kallisto\_plots/SRR950078\_vs\_ERR1831364\_kallisto\_4gb.pdf
