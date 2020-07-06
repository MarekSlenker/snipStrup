#!/bin/bash

# Author: Marek Šlenker
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html


while getopts "r:d:s:v:f:" INITARGS; do
	case "$INITARGS" in
		r) REFSEQ="$OPTARG";;  # reference sequence
		d) DATADIR="$OPTARG";; # path to a directory containing forward and reverse FASTQ files
		v) VCFDIR="$OPTARG";; # output directory for vcf files
		f) SAMPLEPLOIDYFILE="$OPTARG";; # path to a file containing sample names and their ploidy in expected structure
	esac
done




echo "Calling SNPs..."
echo

while read SAMPLE; do 
	read PLOIDY;
	REFSEQBN="$(basename "$REFSEQ")"
	
	echo "** Processing $SAMPLE with ploidy = $PLOIDY"
	echo

	## STEP ONE: Map reads
	echo "Mapping Reads"
	
	bwa mem "$REFSEQ" "$DATADIR"/"$SAMPLE".R1.fq "$DATADIR"/"$SAMPLE".R2.fq -R "@RG\tID:4\tLB:lib1\tPL:ILLUMINA\tSM:$SAMPLE" > "$VCFDIR"/"$SAMPLE".${REFSEQBN%.*}.sam
	samtools sort "$VCFDIR"/"$SAMPLE".${REFSEQBN%.*}.sam -o "$VCFDIR"/"$SAMPLE".${REFSEQBN%.*}.sorted.bam
	
	samtools index "$VCFDIR"/"$SAMPLE".${REFSEQBN%.*}.sorted.bam


	## STEP TWO: Identify variants
	echo "Identifying variants"

	gatk --java-options "-Xmx4g" HaplotypeCaller  \
		-R "$REFSEQ" \
		-I "$VCFDIR"/"$SAMPLE".${REFSEQBN%.*}.sorted.bam \
		-O "$VCFDIR"/"$SAMPLE".${REFSEQBN%.*}.raw.vcf \
		-ploidy "$PLOIDY"
   
	gatk VariantFiltration \
		-R "$REFSEQ" \
		-V "$VCFDIR"/"$SAMPLE".${REFSEQBN%.*}.raw.vcf \
		-O "$VCFDIR"/"$SAMPLE".${REFSEQBN%.*}.filtered.vcf \
		--filter-expression 'QD < 2.0' --filter-name QDfilter \
		--filter-expression 'DP < 8.0' --filter-name DPfilter \
		--filter-expression 'MQ < 40.0' --filter-name MQfilter \
		--filter-expression 'FS > 60.0' --filter-name FSfilter 
	
	## STEP THREE:
	echo cleaning
	rm "$VCFDIR"/*.sam "$VCFDIR"/*.bam "$VCFDIR"/*.idx "$VCFDIR"/*.bai
	
	  
	echo "** $SAMPLE Done"
	echo "******************"

done < "$SAMPLEPLOIDYFILE"


echo
echo "SNPs done!"
echo

exit


