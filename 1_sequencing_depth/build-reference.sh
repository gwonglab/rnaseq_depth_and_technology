#!/bin/bash

#SBATCH -J build_ref	# set a job name
#SBATCH -N 1		# request 1 node for job
#SBATCH -n 1		# request 1 task
#SBATCH -c 1		# request 1 cpu per task
#SBATCH --export=PATH
#SBATCH --mem=10000     # 7670 megabytes were insufficient for some species

# usage: build-reference.sh <genomic-fasta.fna.gz>

# convert Fasta genome and GFF transcriptome (both Gzip compressed) into
# pre-build indicies used by Bowtie2, Tophat2, and Blat 
# assumes input files named with pattern basename.fna.gz and basename.gff.gz
# output files placed in a directory named  REF in the current directory


CASE="$1"	     # assign command-line argument to $CASE
REFDIR="$(pwd)/REF"

# build all the reference files in $REFDIR - create if nonexistant 
BASE=$(basename "${CASE}" .fna.gz)
mkdir -p "${REFDIR}"

# uncompress and change extension to .fa (required by bowtie2)
gzip -d -c "${CASE}" > "${REFDIR}/${BASE}.fa"
gzip -d -c "${CASE%.fna.gz}.gff.gz" > "${REFDIR}/${BASE}.gff"

# build a Bowtie2 index of the genome - needed by Tophat2
bowtie2-build "${REFDIR}/${BASE}.fa" "${REFDIR}/${BASE}"

# build a Tophat index of the transcriptome - saves time
tophat2 --transcriptome-index="${REFDIR}/${BASE}/${BASE}" \
   --GTF="${REFDIR}/${BASE}.gff" "${REFDIR}/${BASE}"

# build a Blat common k-mer database
blat -makeOoc="${REFDIR}/${BASE}-11.ooc" "${REFDIR}/${BASE}.fa" /dev/null /dev/null

