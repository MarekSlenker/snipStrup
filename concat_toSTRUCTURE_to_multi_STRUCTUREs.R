
# Author: Marek Šlenker
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# args: 
# [1] toSTRUCTUREDIR
# [2] 1 or 0, to specify whether to write PopData columns to STRUCTURE file
# [3] a comma-separated string specifying population's membership of samples
# [4] MAXPLOIDY - desired ploidy of structure file
# [5] name prefix of output files     
# [6] number of output files  
# [7] How many random SNPs do we take from one toStructure file? e.g. 1, 2, 3... 
    
## Do not exit on error
options(error=expression(NULL))

args <- commandArgs(TRUE) 
args


toStructureFiles = list.files(path = args[1], pattern = ".toSTRUCTURE", full.names = TRUE, recursive = FALSE)

# Make first columns with sample names and pops
firstFile = read.delim(file = toStructureFiles[1], header = F, sep = " ")
XXsamples = as.character(lapply(firstFile[,1][-1], as.character))

if (args[2] == "1") { # POPDATA = "1"
  pops = strsplit(x = args[3], fixed = T, split = ",")
  popsXploidy=NULL
  for (i in 1:length(pops[[1]])) {
    popsXploidy = c(popsXploidy, rep(pops[[1]][i], as.numeric(args[4])))
  }
  samples = data.frame(XXsamples,  popsXploidy)
  } else{  # POPDATA != "1"
  samples = as.data.frame(XXsamples)
}
samplesPrevious=XXsamples

# prepare STRUCTURE files
for (i in 1:as.numeric(args[6])) {
  assign(paste("structureFile",i,sep = ""), samples)
}



# fill STRUCTURE files with data

for (file in 1:length(toStructureFiles)) {
    
  cat(paste("Processing file",file, "out of", length(toStructureFiles),": ", toStructureFiles[file],"\n"))
  
  nextFile = read.delim(file = toStructureFiles[file], header = F, sep = "")

  snpNames = as.character(lapply(nextFile[1,][- length(nextFile[1,])], as.character))
  snpValues = nextFile[-1,-1]

  samplesCurrent = as.character(lapply(nextFile[,1][-1], as.character))

  if (!all(samplesPrevious==samplesCurrent)) {stop(paste("samples of ",toStructureFiles[file]," are in different order than previous", sep = ""))
  } else { 
      samplesPrevious = samplesCurrent }

  colnames(snpValues) = snpNames
  
  if (length(snpValues) >= as.numeric(args[6])*as.numeric(args[7])) {
    snpsToTake = sample(1:length(snpValues), as.numeric(args[6])*as.numeric(args[7]), replace=FALSE)
  } else {
    snpsToTake = sample(1:length(snpValues), as.numeric(args[6])*as.numeric(args[7]), replace=TRUE)
  }
  
  for (fileN in 1:as.numeric(args[6])) {
    strFile = get(paste("structureFile",fileN,sep = ""))
    snpPos = seq(fileN*as.numeric(args[7]), ((fileN-1)*as.numeric(args[7])+1),-1)
    for (i in 1:length(snpPos)) {
      strFile = cbind(strFile, snpValues[snpsToTake[snpPos[i]]])
    }
    assign(paste("structureFile",fileN,sep = ""),strFile) 
  }
  
}

# remove colnames  
for (fileN in 1:as.numeric(args[6])) {
  strFile = get(paste("structureFile",fileN,sep = ""))
  if (args[2] == "1") { # POPDATA = "1"
    colnames(strFile)[1]=""
    colnames(strFile)[2]=""
  } else{  # POPDATA != "1"
    colnames(strFile)[1]=""
  }
  assign(paste("structureFile",fileN,sep = ""),strFile) 
}

# write files
cat(paste("Writting files...\n"))
for (fileN in 1:as.numeric(args[6])) {
  strFile = get(paste("structureFile",fileN,sep = ""))
  
  write.table(x = strFile, file = paste(args[5], "_", fileN, ".STRUCTURE", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)

}


  
# exit
  
