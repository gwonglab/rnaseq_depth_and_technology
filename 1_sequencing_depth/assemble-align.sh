#!/bin/bash

#SBATCH -J assembly     # set a job name
#SBATCH -N 1            # request 1 node
#SBATCH -n 1            # request 1 task
#SBATCH -c 1            # request 1 cpu per task
#SBATCH --export=PATH
#SBATCH --mem=40000

# usage: sbatch --array=1-6,8-16:2 assembly-full.sh SRR1523365

# Uses subset.py to select a multi-gigabase subset which is then assembled with
# SOAPdenovo-Trans and the scaffolds aligned to the reference with Blat.
# Reference name is looked-up using an internal table.

 
BASE="$1"		      # basename of the dataset
SIZE=${SLURM_ARRAY_TASK_ID})  # size of subset in gigabases 
THREADS=${SLURM_JOB_CPUS_PER_NODE}
REFDIR="$(pwd)/REF"

# create new output directories for the assembly
NUM_GB=$(printf %02dGB ${NUM})
mkdir -p "${BASE}-assembly/${BASE}-${NUM_GB}"
cd       "${BASE}-assembly/${BASE}-${NUM_GB}"


# randomly select a $NUM gigabases of the mapped sequence reads
../../subset.py $((SIZE*1000000000)) "${BASE}-${NUM_GB}" \
    ../../${BASE}-parts/${BASE}-read-stats \
    ../../${BASE}-parts/${BASE}-???/${BASE}-???-mapped-1.fq


# assemble the subset with SOAPdenovo-Trans and GapCloser
# wrap in a command group to allow redirecting STDOUT for logging
{
   # write a config file for the assembly job
   echo -e "max_rd_len=150\n[LIB]\navg_ins=250\nrank=1\nq1=${BASE}-${NUM_GB}-1.fq\nq2=${BASE}-${NUM_GB}-2.fq" > "${BASE}-${NUM_GB}.config"

   # run the assembler programs
   SOAPdenovo-Trans-31mer all -s "${BASE}-${NUM_GB}.config" -p ${THREADS} \
      -r -F -o "${BASE}-${NUM_GB}"
   GapCloser -a "${BASE}-${NUM_GB}.scafSeq" -b "${BASE}-${NUM_GB}.config" \
      -o "${BASE}-${NUM_GB}.GapCloser.fa" -l 150 -p 25 -t ${THREADS}

# end of the command group - write to log files
} > "${BASE}-${NUM_GB}.SOAPdenovo-log"



# lookup the species reference using $BASE
REFNAME=$(grep "^${BASE}" <<-EOF | cut -f 2
	DRR018424	GCF_000001735.3_TAIR10_genomic			Arabidopsis thaliana
	SRR1061361	GCF_000001735.3_TAIR10_genomic			Arabidopsis thaliana
	SRR1178906	GCF_000005425.2_Build_4.0_genomic		Oryza sativa japonica
	SRR1509508	GCF_000001215.4_Release_6_plus_ISO1_MT_genomic	Drosophila melanogaster
	SRR1523365	GCF_000002985.6_WBcel235_genomic		Caenorhabditis elegans
	SRR1732347	GCF_000001635.24_GRCm38.p4_genomic		Mus musculus
	SRR1047863	GCF_000001405.29_GRCh38.p3_genomic		Human
	SRR980471	GCF_000001405.29_GRCh38.p3_genomic		Human
EOF
   )

# align the assembled scaffolds to the reference genome
REFBASE="${REFNAME#GCF_}"
PSL="${BASE}-${NUM_GB}-GCF_${REFBASE%%_*}"
blat "${REFDIR}/${REFNAME}.fa" "${BASE}-${NUM_GB}.GapCloser.fa" \
    -ooc="${REFDIR}/${REFNAME}-11.ooc" -fine "${PSL}-fine.psl"

# # filter the alignments - only high quality alignments - max 1 per scaffold
# # tail skips table header; sort by scaffold name, then by aligned bases
# # awk checks aligned bases against length & scaffold name against prior output 
# tail -n +6 "${PSL}-fine.psl" \
#    | sort -k 10.2,10.2r -k 10.10,10n -k 10.3,10n -k 1,1nr \
#    | awk '($1 >= 0.95 * $11 && $10 != PRIOR) {print; PRIOR = $10}' \
#    > "${PSL}-fine-filter.table"

