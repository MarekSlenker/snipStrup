
# Author: Marek Šlenker
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# args: 
# [1] toSTRUCTUREDIR
# [2] 1 or 0, to specify whether to write PopData columns to STRUCTURE file
# [3] a comma-separated string specifying population's membership of samples
# [4] MAXPLOIDY - desired ploidy of structure file
# [5] a name of output file



## Do not exit on error
options(error=expression(NULL))

args <- commandArgs(TRUE) 
args


toStructureFiles = list.files(path = args[1], pattern = ".toSTRUCTURE", full.names = TRUE, recursive = FALSE)


firstFile = read.delim(file = toStructureFiles[1], header = F, sep = " ")
cat(paste("Processing file 1 out of", length(toStructureFiles),": ", toStructureFiles[1],"\n"))
  
snpNames = as.character(lapply(firstFile[1,][- length(firstFile[1,])], as.character))
snpValues = firstFile[-1,-1]
XXsamples = as.character(lapply(firstFile[,1][-1], as.character))
colnames(snpValues) = snpNames

if (args[2] == "1") { # POPDATA = "1"
  pops = strsplit(x = args[3], fixed = T, split = ",")
  popsXploidy=NULL
  for (i in 1:length(pops[[1]])) {
    popsXploidy = c(popsXploidy, rep(pops[[1]][i], as.numeric(args[4])))
  }

STRUCTUREFile = cbind(data.frame(XXsamples,  popsXploidy), snpValues)
} else{  # POPDATA != "1"
  STRUCTUREFile = cbind(as.data.frame(XXsamples), snpValues)
}

if (length(toStructureFiles)>1){
   samplesPrevious=XXsamples

   for (file in 2:length(toStructureFiles)) {
  
      cat(paste("Processing file",file, "out of", length(toStructureFiles),": ", toStructureFiles[file],"\n"))
  
      nextFile = read.delim(file = toStructureFiles[file], header = F, sep = " ")

      snpNames = as.character(lapply(nextFile[1,][- length(nextFile[1,])], as.character))
      snpValues = nextFile[-1,-1]

      samplesCurrent = as.character(lapply(nextFile[,1][-1], as.character))

      if (!all(samplesPrevious==samplesCurrent)) {
        stop(paste("samples of ",toStructureFiles[file]," are in different order than previous", sep = ""))
      } else {
        samplesPrevious = samplesCurrent
      }

      colnames(snpValues) = snpNames

      STRUCTUREFile = cbind(STRUCTUREFile, snpValues)
    }
}

if (args[2] == "1") { # POPDATA = "1"
	colnames(STRUCTUREFile)[1]=""
	colnames(STRUCTUREFile)[2]=""
} else{  # POPDATA != "1"
	colnames(STRUCTUREFile)[1]=""
}
  

write.table(x = STRUCTUREFile, file = paste(args[5],".STRUCTURE", sep=""), sep = "\t", row.names = F, col.names = T, quote = F)
  
  
# exit
  
