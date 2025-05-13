printprint <- function(i = '') { print(paste0("print('", i, "')"))}

GenotypeFiltering <- function(inGenotypes, minPresent = 75, maxPresent = 10000, createRDS = FALSE, rdsName = "filteredGenotypes.rds") {
  nLoci <- ncol(inGenotypes) - 1
  sumEmpties <- function(x) {return(sum(grepl("^$", x)))}
  
  filteredGenotypes <- inGenotypes[apply(inGenotypes, 1, sumEmpties) <= (nLoci - minPresent), ]
  filteredGenotypes <- filteredGenotypes[apply(filteredGenotypes, 1, sumEmpties) >= (nLoci - maxPresent), ]
  
  if(createRDS == TRUE) {
    saveRDS(filteredGenotypes, file = rdsName)
  }
  return(filteredGenotypes)
}

### THE FOLLOWING FUNCTIONS ARE DEPRECATED ###

CompareBears <- function(bear1Data, bear2Data, nLoci) { # count matching loci between two bears
  simSum <- sum(bear1Data == bear2Data, na.rm = TRUE) # number of shared alleles between the two samples
  return(nLoci - simSum) # return number of non shared alleles between the two samples
}

BearDist <- function(snpFrame) { # function that generates a dist object with how many loci differ and or are missing between individuals
  sampNames <- snpFrame$Sample
  nSamples <- length(sampNames)
  GetDiffs <- function(current, nSamples, snpFrame, nLoci) {
    print(paste0("Calculating sample ", current, " out of ", nSamples))
    remainingIndices <- (1:nSamples)[(current + 1):nSamples]
    sampleDiffs <- apply(snpFrame[remainingIndices, -c(1)], 1, CompareBears, snpFrame[current, -c(1)], nLoci) # uses experimental more efficient comparison
    return(sampleDiffs)
  }
  snpFrame[snpFrame == ''] <- NA
  nLoci <- ncol(snpFrame) - 1
  diffList <- lapply(1:(nSamples - 1), GetDiffs, nSamples, snpFrame, nLoci)
  bearDiffs <- do.call(c, diffList)
  bearMatrix <- matrix(data = 0, nrow = nSamples, ncol = nSamples)
  colnames(bearMatrix) <- sampNames
  rownames(bearMatrix) <- sampNames
  bearDistObject <- as.dist(bearMatrix)
  bearDistObject[] <- unlist(bearDiffs) # fill dist object with distance values
  return(bearDistObject)
}

### DEPRECATED BY PYTHON METHOD:
FastqLoader <- function(fileName, pathToData = filteredFolder) { # load fastq
  loadedFile <- ShortRead::readFastq(dirPath = pathToData, pattern = fileName)
  print(paste("Searching for", fileName))
  print(paste("Found", length(list.files(pathToData, pattern = fileName)), "files."))
  return(c(loadedFile, fileName))
}

PickGenotype <- function(srVector, cutoff) {
  srObject <- srVector[[1]]
  srTable <- ShortRead::tables(srObject)$top
  srTable <- srTable[names(srTable) != ""] # remove empty "SNP"
  if(length(srTable) > 1) {
    topSNP <- srTable[1]
    sndSNP <- srTable[2]
  } else if(length(srTable) == 1) {
    topSNP <- srTable[1]
    sndSNP <- 0 # as to not duplicate number of reads
  } else {
    topSNP <- 0
    sndSNP <- 0
  }
  if((topSNP + sndSNP) >= cutoff[1]) {
    if(((sndSNP / (sndSNP + topSNP)) * 100) >= cutoff[2]) { # if a heterozygote or if there was only one top nucleotide
      # print(paste0(sort(c(names(topSNP), names(sndSNP))), collapse = ""))
      return(c(srVector[[2]], paste0(sort(c(names(topSNP), names(sndSNP))), collapse = "")))
      # return(c(srVector[[2]], paste0(names(topSNP), names(sndSNP))))
    } else { # if a homozygote
      # return(c(srVector[[2]], names(topSNP), names(topSNP)))
      return(c(srVector[[2]], paste0(names(topSNP), names(topSNP))))
    }
  } else {
    # return(c(srVector[[2]], "-", "-"))
    return(c(srVector[[2]], ""))
    # return(c(srVector[[2]], NA))
  }
}

RunLocus <- function(locusName, allFiles, filteredFolder, locusCutoffs = c()) {
  print(paste0("Running locus: ", locusName))
  locusFiles <- allFiles[grepl(locusName, allFiles)]
  loadedFastqFiles <- lapply(locusFiles, FastqLoader, filteredFolder)
  if(locusName %in% names(locusCutoffs)) { # if a specific cutoff is set for this locus
    cutoff = locusCutoffs[locusName][[1]]
  } else { # else use standard settings
    cutoff = c(10, 1) # minimum 10 reads, minimum 10% to be a heterozygote
  }
  genotypeRes <- lapply(loadedFastqFiles, PickGenotype, cutoff)
  rm(loadedFastqFiles) # remove to conserve memory
  gc() # run garbage collection to free memory
  genotypeCols <- do.call(rbind, genotypeRes)
  genotypeCols <- as.data.frame(genotypeCols)
  # colnames(genotypeCols) <- c("Sample", paste0(locusName, "_1"), paste0(locusName, "_2"))
  colnames(genotypeCols) <- c("Sample", locusName)
  genotypeCols$Sample <- gsub(".+_|BMK.+", "", genotypeCols$Sample)
  return(genotypeCols)
}

RunAllLoci <- function(inputFiles, filteredFolder, locusCutoffs = c(), specifiedLoci = c()) {
  lociNames <- unique(gsub("_.+", "", inputFiles)) # remove everything after and including the first underscore
  genotypesList <- lapply(lociNames, RunLocus, inputFiles, filteredFolder, locusCutoffs)
  allGenotypes <- Reduce(function(x, y) merge(x, y, by = "Sample"), genotypesList)
  if(length(specifiedLoci) > 0) {
    allGenotypes <- allGenotypes[, c("Sample", specifiedLoci)]
  }
  return(allGenotypes)
}

GenerateCutoffs <- function(inputFiles, minReads = 10, minPerc = 20) { # min reads to keep locus, min % to keep lesser allele
  locusCutoffs <- c()
  locusCutoffs <- rep(list(c(minReads, minPerc)), length(unique(gsub("_.+", "", inputFiles))))
  names(locusCutoffs) <- unique(gsub("_.+", "", inputFiles))
  return(locusCutoffs)
}

# The above deprecated functions can be used like this if python is a problem:
#filteredFolder <- "./Filtered_data/SNP_filtered/"
#inputFiles <- list.files(filteredFolder, pattern = ".fastq.gz")
#locusCutoffs <- GenerateCutoffs(minReads = 10, minPerc = 20)
#locusCutoffs$`locus name` <- c(10, 10) # add 'locus name' and its designated cutoffs, optional
#locusNames <- readLines("bear_snp_loci.csv") # if loci are to be in a specific order, load this
#allGenotypes <- RunAllLoci(inputFiles, filteredFolder, locusCutoffs, locusNames)