# snipStrup
### SNPs STRUCTURE Parser  
**Suitable for mixed polyploid data sets!**
&nbsp;  
&nbsp;  

The **snipStrup** is currently under development. Meanwhile, we make available shell and R scripts for calling SNPs, parsing vcf files to toSTRUCTURE files, combining them and producing STRUCTURE output files, as done in [Melichárková et al. 2020](https://www.ncbi.nlm.nih.gov/pubmed/______). 


### Software Dependencies
BWA: http://bio-bwa.sourceforge.net/bwa.shtml  
Picard: https://broadinstitute.github.io/picard/  
GATK: https://gatk.broadinstitute.org/hc/en-us  
SAMtools: http://samtools.sourceforge.net/  
R: https://www.r-project.org/  


Be sure that the programs are in `$PATH` and the environment variable `$PICARD` points to picard.jar file.


### Installation
You can download this repository zipped (button on the right-hand side of the screen) or, if you have git installed on your system, clone it with:

```bash
git clone https://github.com/MarekSlenker/snipStrup.git
```


### Pipeline Input 
* **Sequencing reads in FASTQ format** - we recommend to remove PCR duplicates, poor quality reads, trim poor quality base calls, and remove adapter sequences.  
*expected naming convention*:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    <sample_name>.R1.fq.bz2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    <sample_name>.R2.fq.bz2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Cmat_TUL8_2x.R1.fq.bz2`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Cmat_TUL8_2x.R2.fq.bz2`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Cprat_LOCH1_4x.R2.fq.bz2`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Cprat_LOCH1_4x.R2.fq.bz2`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Cprat_MOJ2_6x.R2.fq.bz2`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Cprat_MOJ2_6x.R2.fq.bz2`  
      

* **reference sequences** - we will align reads to these sequences (keep them in a single file) and look for differences between reads and these reference sequences.  
*expected naming convention*:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    <file_with_ref_sequences>.fasta  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Camara_reference.fasta`  

* **SAMPLEPLOIDYFILE** - a file containing <sample_name> and <sample_ploidy> pairs (in separate lines).  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Cmat_TUL8_2x`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `2`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Cprat_LOCH1_4x`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `4`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Cprat_MOJ2_6x`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `6`  

* **GCPFILE** - Gene Concatenation Pattern File - **not required**. The pipeline produce one `toSTRUCTURE` file for each sequence. If you have multiple exon sequences (multiple `toSTRUCTURE` files) for a single gene, you can concate them to single `toSTRUCTURE` file representing one gene.
File contains <patterns> in separate lines. The `toSTRUCTURE` files containing same pattern in the name will be concatenated.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Assembly_000000003997`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Assembly_000000014029`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Assembly_000000033889`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Assembly_000000055321`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Assembly_000000057671`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    `Assembly_000000062661`  

&nbsp;  

### Pipeline Output
This pipeline will map the reads to the reference sequences, call and filter detected variants, and convert vcf files to toSTRUCTURE files (one per sequence). 
If you have more exon sequences for one gene (as in test data), toSTRUCTURE files of these sequences should be concatenated. We expect that many SNPs will be found on the same reference sequence, thus merging all toSTRUCTURE files to single STRUCTURE input file can violate the assumption of no linkage between sites. Thus, the preferred approach is to produce a reasonable amount of datasets by drawing a random X% SNP sites from each gene (toSTRUCTURE file). These datasets should be analysed separately and the results combined using the program CLUMPP (https://rosenberglab.stanford.edu/clumpp.html). If there is any reason, toSTRUCTURE files can be merged to a single STRUCTURE input file by provided alternative script. The STRUCTURE software requires that each individual takes the same number of rows that is equal to the highest ploidy in the data set. Samples of a lower ploidy than given in the `PLOIDY` parameter (for mixed-ploidy data sets) must have a missing data symbol inserted to fill in the extra rows. This pipeline will end up with ready-to-run STRUCTURE input files. 


## Running the pipeline
This pipeline consists of five parts, each is done by a separate script. The `Main.sh` is a wrapper script connecting four parts to one pipeline. Please modify it according to your needs. If you have access to HPC cluster, we highly recommend splitting the first step (variant calling) to separate tasks - one task for one reference sequence. The implementation in `Main.sh` reaches the end in a reasonable time only with a small dataset and is mainly for demonstration purposes. 

Scripts depend on the following list of variables. You can set them in the `Main.sh` script. Please ensure they are set up correctly.  
* `DATADIR` = path to a directory containing forward and reverse FASTQ files (assuming paired-end reads ) for each individual. For expected naming convention see above.  
* `REFSEQ` = reference fasta sequences. For expected naming convention see above.  
* `SAMPLEPLOIDYFILE`= path to a file containing sample names and their ploidy in expected structure (see above).  
* `GCPFILE`= path to a file containing patterns for concatenation (see above). Not required.  
* `SCRIPTDIR`= path to this directory, containing running scripts.  
&nbsp;
* `VCFDIR` = output directory for vcf files.  
* `toSTRUCTUREDIR` = output directory for toSTRUCTURE files.
* `STRUCTUREDIR` = output directory for STRUCTURE files.  
&nbsp;
* `STRUCTURE_OUTFILES` = number of random datasets to produce - STRUCTURE input files  
* `STRUCTURE_OUTFILE_PREFIX` = a prefix to be added to STRUCTURE output file names.  
* `STRUCTURE_SNPS_TO_OUTFILE` = How many random SNPs should we take from one toStructure file? Percent or integer. Depends on the script.
* `STRUCTURE_POPDATA` = `1` or `0`, to specify whether to write PopData columns to STRUCTURE file. For more details see STRUCTURE manual.  
* `STRUCTURE_POPS` = a comma-separated string specifying population's membership of samples. Eg., "1,2,2,2,3". Not required, but cannot be empty! Used only when `STRUCTURE_POPDATA="1"`.  




As next, we recommended checking default values of hard-filtering criteria in `callSNPs.sh` script. Filtration is done by `VariantFiltration`, filters are applied through `--filter-expression` argument. See GATK documentation for more details.

Now you can run the pipeline typing
```bash
./Main.sh 
```

The `Main.sh` script will do the following:
1. Run `callSNPs.sh` script, which will do the following:
      - Map paired-end reads to the reference sequences using `bwa mem`.
      - Call SNPs using `gatk HaplotypeCaller`.
      - Filter variant calls based on certain criteria using `gatk VariantFiltration`.  
      
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; In the end, you will have vcf files storing variants of each sample.

2. Run `returnHighestPloidyOfSamples.sh` to find the highest ploidy in the dataset, i.e., how many rows per individual to write to STRUCTURE file.

3. Run `VCF_to_toSTRUCTURE.R` to convert sample's vcf files to \"toSTRUCTURE\" file (separate toSTRUCTURE file for each reference sequence will be generated), remove filtered variants, and remove indels as these variants are not supported in STRUCTURE. Please read R log file, as it may contain warnings. 

4. Run `concat_toSTRUCTURE_by_pattern_to_toSTRUCTURE.R` to concat exon's toSTRUCTURE files to gene's toSTRUCTURE fiels. Not required.

5. Run `concat_toSTRUCTURE_to_multi_STRUCTUREs_by_PERCENTS.R` to produce defined numbers of STRUCTURE files, each containing random N% SNP site from each gene. This option ensures taking same portion of SNPs from diversely variable genes and on the other hand eliminates genes with low level of variation. (If you want to take 5%, genes with less than 11 variants will be discarded, because 10 is after rounding less than one).  
Alternatively, to take same amount of SNPs from each toSTRUCTURE file, run `concat_toSTRUCTURE_to_multi_STRUCTUREs.R`, or `concat_toSTRUCTURE_to_single_STRUCTURE.R` to merge all toSTRUCTURE files to single STRUCTURE file.



### Equivocal variant filtration
Variant filtration is performed for each vcf file, and you can edit hard-filtering criteria in `callSNPs.sh` script. However, consider the scenario, in which the same variant is kept in one but filtered out in another sample (regardless of the reason). Because it is more acceptable to remove than to increase arbitrarily the variability, variants that did not pass filtering in some samples are removed from the whole dataset.


### Test Dataset
We provide a test dataset that is a subset of real Illumina data from an enriched library. To try the functionality of this package, set these parameters in `Main.sh` script and run `./Main.sh` from your command line.

```
DATADIR="test_data"
REFSEQ="$DATADIR"/Cardamine_targetSequences.fasta 
SAMPLEPLOIDYFILE="$DATADIR"/samples_ploidy_list.txt 
GCPFILE="$DATADIR"/gene_concatenation_pattern.txt
SCRIPTDIR="."

VCFDIR="vcf_files"             # output directory for vcf sequences
toSTRUCTUREDIR="toSTRUCTURE"   # output directory for toSTRUCTURE files
CONCATENATED="toSTRUCTURE_concatenated" # output directory for concatenated toSTRUCTURE files
STRUCTUREDIR="STRUCTURE_files" # output directory for STRUCTURE files

STRUCTURE_OUTFILES="100"              # number of random datasets to produce - STRUCTURE input files  
STRUCTURE_OUTFILE_PREFIX="testData"
STRUCTURE_SNPS_TO_OUTFILE="0.05"      # How many random SNPs do we take from one toStructure file? percent or integer, depends on script
STRUCTURE_POPDATA="1"              # 1/0  - input file contains a population identifier
STRUCTURE_POPS="1,2,2,2,3"  
```

Parameters in file STRUCTURE's *mainparams* file:
```
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line

#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   1     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier

#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data before the genotype data start.

#define MARKERNAMES      1  // (B) data file contains row of marker names
#define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs) and a row to indicate which alleles are recessive
#define MAPDISTANCES     0  // (B) data file contains row of map distances between loci
```

## Other questions not covered here and reporting problems
If you have a question or you encounter a problem, please see [issues](https://github.com/MarekSlenker/Hyb-Seq_SNPs_STRUCTURE/issues) and feel free to ask any question. The authors will do their best to help you.
