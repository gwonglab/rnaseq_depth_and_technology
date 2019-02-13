#!/bin/bash

#SBATCH -J map-reads	# set a job name
#SBATCH -N 1            # request 1 node
#SBATCH -n 1		# request 1 task
#SBATCH -c 1		# request 1 cpu per task
#SBATCH --mem=5700

# usage: sbatch --array=1-324 ./map-reads.sh SRR1523365.sra
# 
# Uses fastq-dump (SRA Toolkit) to extract one million reads from a SRA file
# before trimming primer/adapter sequence and mapping to a reference
# transcriptome/genome using Tophat2.
# The primer/adapters and transcriptomes are in look-up tables in this script
# based on the SRA filename.


SRA_FILE=$(realpath "$1") 
BASE="$(basename "$1" .sra)"
PART_NUM=$SLURM_ARRAY_TASK_ID
THREADS=$SLURM_JOB_CPUS_PER_NODE
REFDIR="$(pwd)/REF"


# create subdirectories for this output
export LC_ALL=C	  # speed up sorts, etc. by simplifying ordering rules
ZERO_PART_NUM=$(printf %03d $PART_NUM)  # zero-pad part no. to three digits
CASE="${BASE}-${ZERO_PART_NUM}"
mkdir -p "${BASE}-parts/${CASE}"
cd       "${BASE}-parts/${CASE}"

# extract a million spots from the SRA file in Fastq format
END=$(($PART_NUM * 1000000))
START=$((END - 999999))
fastq-dump --minSpotId $START --maxSpotId $END \
   --defline-seq '@$si/$ri' --defline-qual + --split-files "$FILE"


# lookup the adapter/primer sequences for the 1st and 2nd ends
AD1=$(grep "^${BASE}" <<-EOF | cut -f 2
	DRR018424	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
	SRR1061361	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
	SRR1523365	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
	SRR1732347	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
	SRR980471	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
	SRR1178906	AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
	SRR1047863	AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
	SRR1509508	AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
EOF
   )
AD2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT



# trim and primer/adapter sequence - no base quality-score based trimming
~/mapping/trim_galore --phred33 --paired --quality 0 \
   -a $AD1 -a2 $AD2 ${BASE}_?.fastq

# # alternative trims 5'-most bases from reads - used for SRR1523365
# ~/mapping/trim_galore --phred33 --paired --quality 0 \
#   --clip_R1 5 --clip_R2 5 -a $AD1 -a2 $AD2 ${BASE}_?.fastq



# use the BASE (from SRA filename) to choose the species reference
REFNAME=$(grep "^${BASE}" <<-EOF | cut -f 2
	DRR018424	GCF_000001735.3_TAIR10_genomic			Arabidopsis thaliana
	SRR1047863	GCF_000001405.29_GRCh38.p3_genomic		Human
	SRR1061361	GCF_000001735.3_TAIR10_genomic			Arabidopsis thaliana
	SRR1178906	GCF_000005425.2_Build_4.0_genomic		Oryza sativa japonica
	SRR1509508	GCF_000001215.4_Release_6_plus_ISO1_MT_genomic	Drosophila melanogaster
	SRR1523365	GCF_000002985.6_WBcel235_genomic		Caenorhabditis elegans
	SRR1732347	GCF_000001635.24_GRCm38.p4_genomic		Mus musculus
	SRR980471	GCF_000001405.29_GRCh38.p3_genomic		Human
EOF
   )



# map the reads with Tophat2
tophat2 --num-threads $THREADS --b2-very-sensitive -o "${CASE}-tophat2" \
    --transcriptome-index="${REFDIR}/${REFNAME}/${REFNAME}" \
    "${REFDIR}/${REFNAME}" "${BASE}"_?_val_?.fq


# extract the paired-end alignments as FastQ files
# flag bitfields + 16 = 00011X00X1 or 00101X00X1 (X bits are don't care)
# note that tophat2 does not reliably set bit 0x2 (segments properly aligned)

samtools view "${CASE}-tophat2/accepted_hits.bam" \
   | gawk '(and($2 + 16, 1005) == 97)  {print "@" $1 "/1\t" $10 "\t+\t" $11}' \
   | sort --buffer-size=1G -k 1,1 | tr '\011' '\012' > "${CASE}-mapped-1.fq"
samtools view "${CASE}-tophat2/accepted_hits.bam" \
   | gawk '(and($2 + 16, 1005) == 161) {print "@" $1 "/2\t" $10 "\t+\t" $11}' \
   | sort --buffer-size=1G -k 1,1 | tr '\011' '\012' > "${CASE}-mapped-2.fq"

# count reads and bases in each of the outputs and record in read-stats file
for FILE in "${CASE}-mapped-1.fq" "${CASE}-mapped-2.fq"; do
   F=$(basename "${FILE}" .fq); F=${F/-mapped-/ }; F=${F/-/ };
   cat "$FILE" | paste - - - - \
       | awk "{B+=length(\$2)}; END {print \"$F\", NR, B}"
        >> ../${BASE}-read-stats
done

exit

