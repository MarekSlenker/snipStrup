
# Author: Marek Šlenker
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# args: 
# [1] concat pattern list
# [2] toSTRUCTURE directory
# [3] dir for concatened files

options(error=expression(NULL))
options(warn=1)

args <- commandArgs(TRUE) 
args


if (file.exists(args[1])) {
  toConcatPattern = readLines(con = args[1])
} else {
  stop(paste("cannot open file '",args[1],"'",sep = ""))
}

if (length(toConcatPattern) == 0) {
  stop("file is empty")
}

for (filePttern in toConcatPattern) {
  filesWithPattern = list.files(path = args[2], pattern = filePttern, full.names = TRUE, recursive = FALSE)
  cat(paste("Processing pattern",which(filePttern == toConcatPattern), "out of", length(toConcatPattern),": ", filePttern,"\n"))
  
  if (length(filesWithPattern) == 0) { 
    warning(paste("Skipping pattern ",filePttern,". No files found.", sep="") )  
    next
  }
  
  if (length(filesWithPattern) == 1) { # only one file
    
    file = read.delim(file = filesWithPattern, header = F, sep = " ")
    
    fileSamples = as.character(file[,1])[-1]
    fileSnpNames = as.character(lapply(file[1,][- length(file[1,])], as.character))
    filesnpValues = as.data.frame(file[-1,-1])
    
    colnames(filesnpValues) = fileSnpNames
    filesnpValues=t(filesnpValues) # musim to teraz drbnut naopak, lebo by to malo problem s rovnakymi rownames
    colnames(filesnpValues) = fileSamples # toto je OK, lebo to je drbnute naopak.. :)
    
    # wite concatened file
    write.table(x = t(filesnpValues), file =  paste(args[3],"/", filePttern, ".toSTRUCTURE", sep = "") , row.names = TRUE, col.names = TRUE)
    
    rm(file, fileSamples, fileSnpNames, filesnpValues)
  } else { # more then one file
    
    firstFile = read.delim(file = filesWithPattern[1], header = F, sep = " ")
    
    samples = as.character(firstFile[,1])[-1]
    firstFileSamples = samples
    
    snpNames = as.character(lapply(firstFile[1,][- length(firstFile[1,])], as.character))
    snpValues = as.data.frame(firstFile[-1,-1])
    
    for (fiWiPatt in 2:length(filesWithPattern)) {
      
      nexttFile = read.delim(file = filesWithPattern[fiWiPatt], header = F, sep = " ")
      
      nextFilesamples = as.character(nexttFile[,1])[-1]
      
      if (!all(firstFileSamples == nextFilesamples)) {
        stop(paste("samples of file",filesWithPattern[fiWiPatt]," are in different order.", sep = ""))
      }
      
      # nepredpokladam ze dlzka je 0, taketo su vyhadzovane o skript skôr
      nexttFileSnpNames = as.character(lapply(nexttFile[1,][- length(nexttFile[1,])], as.character))
      nexttFileSnpValues = as.data.frame(nexttFile[-1,-1])
      
      snpNames=c(snpNames, nexttFileSnpNames)
      snpValues=cbind(snpValues, nexttFileSnpValues)
    }
    
    # make toSTRUCTURE table
    colnames(snpValues) = snpNames
    snpValues=t(snpValues) # musim to teraz drbnut naopak, lebo by to malo problem s rovnakymi rownames
    colnames(snpValues) = samples # toto je OK, lebo to je drbnute naopak.. :)
      
    # wite concatened file
    write.table(x = t(snpValues), file =  paste(args[3],"/", filePttern, ".toSTRUCTURE", sep = "") , row.names = TRUE, col.names = TRUE)
    
    # rm
    rm(firstFile, samples, firstFileSamples, snpNames, snpValues, fiWiPatt, nexttFile, nextFilesamples, nexttFileSnpNames, nexttFileSnpValues)
  }
  
  rm(filesWithPattern)
}





