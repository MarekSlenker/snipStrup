#!/bin/bash

# Author: Marek Šlenker
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# A wrapper script connecting remaining four scripts to one pipeline. Please modify it according to your needs.



# Set initial variables
DATADIR="test_data"   # Fastq data to process
REFSEQ="$DATADIR"/Cardamine_targetSequences.fasta     # reference fasta sequences
SAMPLEPLOIDYFILE="$DATADIR"/samples_ploidy_list.txt     # file storing sample names and their ploidy, for examples see "test_dataset"
GCPFILE="$DATADIR"/gene_concatenation_pattern.txt    # path to a file containing patterns for concatenation
SCRIPTDIR="."                # directory with scripts

# Directories for the partial results
VCFDIR="vcf_files"             # output directory for vcf sequences
toSTRUCTUREDIR="toSTRUCTURE"   # output directory for toSTRUCTURE files
CONCATENATED="toSTRUCTURE_concatenated" # output directory for concatenated toSTRUCTURE files
STRUCTUREDIR="STRUCTURE_files" # output directory for STRUCTURE files


# STRUCTURE parameters
STRUCTURE_OUTFILES="100"              # number of random datasets to produce - STRUCTURE input files  
STRUCTURE_OUTFILE_PREFIX="testData"
STRUCTURE_SNPS_TO_OUTFILE="1"      # How many random SNPs do we take from one toStructure file? percent or integer, depends on script
STRUCTURE_POPDATA="1"              # 1/0  - input file contains a population identifier
STRUCTURE_POPS="1,2,2,2,3"         # values for POPDATA





#######################################################################################################


mkdir "$VCFDIR" 
mkdir "$toSTRUCTUREDIR"
mkdir "$STRUCTUREDIR"
mkdir "$CONCATENATED"

# Decompress FASTQ files
bunzip2 -v "$DATADIR"/*.bz2


## STEP ONE: Process all reference files, call SNPs on the basis of these files. 

# Create a sequence dictionary and indexes
java -jar $PICARD CreateSequenceDictionary R="$REFSEQ"
bwa index "$REFSEQ"
samtools faidx "$REFSEQ"
	
"$SCRIPTDIR"/callSNPs.sh -r "$REFSEQ" -d "$DATADIR" -v "$VCFDIR" -f "$SAMPLEPLOIDYFILE"



## STEP TWO: we need to know the highest ploidy of samples, to reserve these rows in STRUCTURE file
MAXPLOIDY=$( $SCRIPTDIR/returnHighestPloidyOfSamples.sh $SAMPLEPLOIDYFILE )
echo "max ploidy of samples: $MAXPLOIDY"



## STEP THREE: we convert sample's vcf files for each reference sequence to \"toSTRUCTURE\" file
echo "converting sample's vcf files for each reference sequence to \"toSTRUCTURE\" file"

R CMD BATCH --no-save --no-restore "--args 
	$REFSEQ
	$SAMPLEPLOIDYFILE
	$MAXPLOIDY
	$VCFDIR
	$toSTRUCTUREDIR" \
	VCF_to_toSTRUCTURE.R VCF_to_toSTRUCTURE.R.log
	
echo "toSTRUCTURE files were writed to $toSTRUCTUREDIR"
echo "WARNING: Previous script skipped $(grep -c "WARNING: Skipping chromosome " VCF_to_toSTRUCTURE.R.log) sequences. Read VCF_to_toSTRUCTURE.R.log for details."
echo


## STEP FOUR: concat exon's toSTRUCTURE files by pattern to gene's toSTRUCTURE files.
echo "Concatening exon's \"toSTRUCTURE\" files to gene's \"toSTRUCTURE\" files"

R CMD BATCH --no-save --no-restore "--args 
	$GCPFILE
	$toSTRUCTUREDIR
	$CONCATENATED" \
	concat_toSTRUCTURE_by_pattern_to_toSTRUCTURE.R concat_toSTRUCTURE_by_pattern_to_toSTRUCTURE.R.log

echo "Concatened toSTRUCTURE files were writed to $CONCATENATED"
echo "WARNING: Previous script skipped $(grep -c "WARNING: Skipping " concat_toSTRUCTURE_by_pattern_to_toSTRUCTURE.R.log) patterns. Read concat_toSTRUCTURE_by_pattern_to_toSTRUCTURE.R.log for details."
echo



## STEP FIVE: merge all toSTRUCTURE files to final STRUCTURE files


# default approach: take a random n (STRUCTURE_SNPS_TO_OUTFILE) of SNP sites from each of the toSTRUCTURE files
R CMD BATCH --no-save --no-restore "--args 
	$toSTRUCTUREDIR
	$STRUCTURE_POPDATA
	$STRUCTURE_POPS
	$MAXPLOIDY
	$STRUCTUREDIR/$STRUCTURE_OUTFILE_PREFIX
	$STRUCTURE_OUTFILES
	$STRUCTURE_SNPS_TO_OUTFILE" \
	concat_toSTRUCTURE_to_multi_STRUCTUREs.R concat_toSTRUCTURE_to_multi_STRUCTUREs.R.log



echo "STRUCTURE files were writed to $STRUCTUREDIR"


exit


# if you want to merge all SNPs to a single output file, uncomment this alternative script
R CMD BATCH --no-save --no-restore "--args 
	$toSTRUCTUREDIR
	$STRUCTURE_POPDATA
	$STRUCTURE_POPS
	$MAXPLOIDY
	$STRUCTUREDIR/$STRUCTURE_OUTFILE_PREFIX" \
	concat_toSTRUCTURE_to_single_STRUCTURE.R concat_toSTRUCTURE_to_single_STRUCTURE.R.log

# if you want to take a random n % of SNP sites (same portion) from each of the toSTRUCTURE files
R CMD BATCH --no-save --no-restore "--args 
	$CONCATENATED
	$STRUCTURE_POPDATA
	$STRUCTURE_POPS
	$MAXPLOIDY
	$STRUCTUREDIR/$STRUCTURE_OUTFILE_PREFIX
	$STRUCTURE_OUTFILES
	$STRUCTURE_SNPS_TO_OUTFILE" \
	concat_toSTRUCTURE_to_multi_STRUCTUREs_by_PERCENTS.R concat_toSTRUCTURE_to_multi_STRUCTUREs_by_PERCENTS.R.log

