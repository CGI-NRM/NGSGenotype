printprint <- function(i = '') { print(paste0("print('", i, "')"))}

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

CompareBears <- function(bear1, bear2, snpFrame) { # count matching loci between two bears
  CheckDiff <- function(allelePair) {
    x <- allelePair[1]
    y <- allelePair[2]
    if((x == y) && (x != "") && (y != "")) { # if alleles are nonempty and the same
      return(1) # count one similarity
    } else {
      return(0) # count zero similarities
    }
  }
  similarities <- lapply(rbind(snpFrame[bear1, ], snpFrame[bear2, ])[, -c(1)], CheckDiff)
  simSum <- sum(unlist(similarities)) # number of shared alleles between the two samples
  return((ncol(snpFrame) - 1) - simSum) # return number of non shared alleles between the two samples
}

BearDist <- function(snpFrame) { # function that generates a dist object with how many loci differ and or are missing between individuals
  sampNames <- snpFrame$Sample
  nSamples <- length(sampNames)
  bearDiffs <- c()
  for(current in 1:(nSamples - 1)) {
    remainingIndices <- (1:nSamples)[(current + 1):nSamples]
    bearDiffs <- c(bearDiffs, lapply(remainingIndices, CompareBears, current, snpFrame))
    # print(paste0(current, ":", remainingIndices))
  }
  bearMatrix <- matrix(data = 0, nrow = nSamples, ncol = nSamples)
  colnames(bearMatrix) <- sampNames
  rownames(bearMatrix) <- sampNames
  bearDistObject <- as.dist(bearMatrix)
  bearDistObject[] <- unlist(bearDiffs) # fill dist object with distance values
  return(bearDistObject)
}