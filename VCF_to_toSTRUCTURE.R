
# Author: Marek Šlenker
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# args:
# [1] path to reference sequence
# [2] Sample ploidy file
# [3] MAXPLOIDY - desired ploidy of structure file
# [4] VCFDIR - dir with vcf files, result of callSNPs.sh script
# [5] output dir


## Do not exit on error
options(error=expression(NULL))


args <- commandArgs(TRUE) # FASTA seq,  VCF postfix,  output, ploidy
args

library(package="memuse", lib.loc="./rpackages")
library(package="ape", lib.loc="./rpackages")
library(package="vcfR", lib.loc="./rpackages")


samples_ploidy_list = read.delim(args[2], header = F)
# split samples_ploidy_list to sample and ploidy tables
samples_ploidy = samples_ploidy_list[as.numeric(row.names(samples_ploidy_list))%%2 == 0, ]
samples_names = samples_ploidy_list[as.numeric(row.names(samples_ploidy_list))%%2 == 1, ]
rm(samples_ploidy_list)


# Load REF sequences
# REF seq -> table
refSequences = read.dna(file = args[1], format="fasta")

if (is.list(refSequences)){
  .rsl = length(refSequences)

  for (r in 1: .rsl) {
    cat(paste("processing sequence ",r,"of the",.rsl,"\n"))
    refsq = as.data.frame(as.character(refSequences[r]))
    assign(gsub("[[:space:]]", "", names(refSequences[r])), refsq)
    rm(refsq)
  }
  rm(.rsl)
} else {
  seqs = (attr(refSequences,"dimnames")[[1]])
  .rsl = length(seqs)

  for (r in 1: .rsl) {
    cat(paste("processing sequence ",r,"of the",.rsl,"\n"))
    refsq = as.character(refSequences[r,])
    assign(gsub("[[:space:]]", "", seqs[r]), refsq[1,])
    rm(refsq)
  }
}

# Load vcf files
chromsList=NULL
for (sample in samples_names) {
  vcfFile=read.vcfR(file=paste(args[4],"/", list.files(path=args[4],pattern=paste(sample,"[0-9A-Za-z._-]*filtered.vcf",sep="")),sep=""))
  chromsList = c(chromsList, getCHROM(vcfFile))
  
  assign(paste(sample,"_allels",sep = ""), extract.gt(vcfFile, element = "GT", as.numeric = FALSE, return.alleles = T))
  assign(paste(sample,"_refAllels",sep = ""), getREF(vcfFile))
  assign(paste(sample,"_filter",sep = ""), getFILTER(vcfFile))
  assign(paste(sample,"_allelePoss",sep = ""), getPOS(vcfFile))
  assign(paste(sample,"_chromPossitions",sep = ""), getCHROM(vcfFile))

  rm(vcfFile)
  cat(paste("*** Loaded VCF file for sample:",sample,"\n"))
}
chromsList = unique(chromsList)

skipChromosome=FALSE
for (chrom in chromsList) {
  cat(paste("Processing chromosome",which(chrom == chromsList), "out of", length(chromsList),": ", chrom,"\n"))
  
  refsq = get(chrom)
  
  # variant calling table
  vcfTable = NULL
  snpPositins = NULL
  possitionsToRemove = NULL
  
  for (sample in samples_names) {
    
    # ploidy
    pld = as.numeric(as.character(samples_ploidy))[which(samples_names == sample)]
    
    chromPossitions=get(paste(sample,"_chromPossitions",sep = ""))
    
    allels = get(paste(sample,"_allels",sep = ""))[which(chromPossitions == chrom)]
    refAllels = get(paste(sample,"_refAllels",sep = ""))[which(chromPossitions == chrom)]
    filter = get(paste(sample,"_filter",sep = ""))[which(chromPossitions == chrom)]
    allelePoss = get(paste(sample,"_allelePoss",sep = ""))[which(chromPossitions == chrom)]
    
    # SNPs table
    snpTable = as.matrix(refsq)
    colnames(snpTable) = sample
    
    if (length(allels) == 0) {
      cat(paste("WARNING: sample", sample,"has 0 alleles.\n"))
      skipChromosome=TRUE
      next
    }
    
    # find & remember positions with: deletions, insertions and bases not passed filtering
    possitionsToRemove = c(possitionsToRemove, allelePoss[which(nchar(refAllels) > 1)]) # deletion in REF
    possitionsToRemove = c(possitionsToRemove, allelePoss[which(nchar(allels)>(pld*2-1))]) # insertion to ALT
    possitionsToRemove = c(possitionsToRemove, allelePoss[which(! grepl("^[CTAGctag/]+$",allels))]) # weard character
    possitionsToRemove = c(possitionsToRemove, allelePoss[which(filter != "PASS")]) # filtering
    
    # remember SNPs positions
    snpPositins = c(snpPositins, allelePoss)
    
    # move SNPs from SNP table to vcfTable
    snpTable[allelePoss] = allels
    
    vcfTable = cbind(vcfTable,snpTable)
    
    rm(pld, allels, refAllels, filter, allelePoss, snpTable)
  }
  rm(refsq, sample)
  
  # output is vcfTable
  
  if (skipChromosome) {
    skipChromosome=FALSE
    cat(paste("WARNING: Skipping chromosome", chrom,"\n\n"))
    next
  }
  
  uniquesSnpPositins = as.numeric(unique(snpPositins))
  uniquesSnpPositins=sort(uniquesSnpPositins)
  
  uniquePossitionsToRemove = as.numeric(unique(possitionsToRemove))
  uniquePossitionsToRemove=sort(uniquePossitionsToRemove)
  
  snpFreqtable = as.data.frame(table(snpPositins))
  onlyInOneSamplePositions = as.numeric(as.character(snpFreqtable$snpPositins[which(snpFreqtable$Freq == 1)]))
  
  
  # mark positoins to keep = uniquesSnpPositins - uniqueFilteredPositions
  vcfTable = cbind(vcfTable, keep=FALSE)
  vcfTable = as.data.frame(vcfTable)
  
  levels(vcfTable$keep) = c(FALSE,TRUE)
  
  vcfTable$keep[uniquesSnpPositins] = "TRUE"
  vcfTable$keep[uniquePossitionsToRemove] = "FALSE"
  vcfTable$keep[onlyInOneSamplePositions] = "FALSE"
  
  
  SNPs = vcfTable[vcfTable$keep == TRUE,]
  SNPs = SNPs[,- which(colnames(SNPs) == "keep")]
  SNPs = drop(SNPs)
  SNPs=t(SNPs)
  # SNPs = SNPs[ order(row.names(SNPs)), ]   BUG
  
  #
  # Convert to structure
  #
  # table is ready, convert it to structure format
  
  
  # riadkov tolko kolko snipov, stlpcov tolko kolko vzoriek * vysledna ploidia. (tabulka je opat transfomovana)
  STRUCTploidy = as.numeric(args[3])
  STRUCTUREtable=matrix(data = NA, nrow = dim(SNPs)[2], ncol = dim(SNPs)[1]*STRUCTploidy)
  colnames(STRUCTUREtable) = rep(row.names(SNPs), each=STRUCTploidy)
  row.names(STRUCTUREtable) = paste("SNP",colnames(SNPs), sep = "")
  
  #
  
  for (position in 1:length(colnames(SNPs))) {
    
    for (sample in 1:length(row.names(SNPs))) {
      
      snps = strsplit(x = SNPs[sample,position], fixed = T, split = "/")
      samplePloidy = as.numeric(as.character(samples_ploidy))[sample]
      
      for (i in 1:STRUCTploidy) {
        
        if (i>length(snps[[1]])) {
          
          if (i>samplePloidy) {
            STRUCTUREtable[position,STRUCTploidy*sample-STRUCTploidy+i] = "-9"
          } else {  # there is only one nucleotide
            STRUCTUREtable[position,STRUCTploidy*sample-STRUCTploidy+i] = snps[[1]]
          }
          
        } else{
          STRUCTUREtable[position,STRUCTploidy*sample-STRUCTploidy+i] = snps[[1]][i]
        }
      }
      
    }
  }
  
  #(A=1, T=2, G=3, C=4)!
  STRUCTUREtable[STRUCTUREtable=="A" | STRUCTUREtable=="a"] = 1
  STRUCTUREtable[STRUCTUREtable=="T" | STRUCTUREtable=="t"] = 2
  STRUCTUREtable[STRUCTUREtable=="G" | STRUCTUREtable=="g"] = 3
  STRUCTUREtable[STRUCTUREtable=="C" | STRUCTUREtable=="c"] = 4
  
  write.table(x = t(STRUCTUREtable),file = paste(args[5], "/", chrom, ".toSTRUCTURE", sep = "") , row.names = T, col.names = T)
  cat("\n")
}

# exit
