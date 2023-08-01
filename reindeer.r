library(ShortRead)

ParsePrimers <- function(fastaFile) {
  fastaString <- paste(readLines(fastaFile), collapse = "\n")
  fastaList <- strsplit(fastaString, split = ">")
  fastaList <- fastaList[[1]][fastaList[[1]] != ""]
  newLineSplit <- function(inputLine) {
    sequenceList <- strsplit(inputLine, split = "\n")[[1]]
    nucSequence <- sequenceList[2:length(sequenceList)]
    nucSequence <- paste(nucSequence, collapse = "")
    sequenceList <- c(sequenceList[1], strsplit(nucSequence, "\\.\\.\\."))
    sequenceDF <- data.frame(Primer = c(sequenceList[[2]][1], sequenceList[[2]][2]))
    colnames(sequenceDF) <- sequenceList[[1]]
    return(sequenceDF)
  }
  fastaDFList <- lapply(fastaList, newLineSplit)
  fastaDF <- do.call("cbind", fastaDFList)
  return(fastaDF)
}

primerDF <- ParsePrimers("./primers.fa")
inputFiles <- list.files("./Filtered_data/", pattern = "*.fastq.gz")

NameShortener <- function(name) {
  # name <- substr(name, 1, 20)
  name <- gsub("_[^_]+$", "", name)
  return(name)
}

BarCol <- function(Allelelength, highlight) {
  c1 <- rep(wesanderson::wes_palette("Chevalier1")[3], length(Allelelength))
  c1[which(names(Allelelength) %in% highlight)] <- wesanderson::wes_palette("Chevalier1")[2]
  c1
}

inputFiles <- unlist(lapply(inputFiles, NameShortener))
inputFiles <- unique(inputFiles)

FastqLoader <- function(fileName) {
  namePattern <- paste0(fileName, "*")
  loadedFile <- readFastq(dirPath = "./Filtered_data/", pattern = namePattern)
  print(paste("Searching for", namePattern))
  print(list.files("./Filtered_data/", pattern = namePattern))
  return(c(loadedFile, fileName))
}

loadedFastqFiles <- lapply(inputFiles, FastqLoader)

FragLength2 <- function(sample, primerF = "AAAAA", primerR = "AAAAA", range = 0:500) {
  lengthF <- nchar(primerF)
  lengthR <- nchar(primerR)
  fragl <- ShortRead::width(ShortRead::sread(sample))
  lengthTot <- lengthF + lengthR + fragl
  lengthTot <- lengthTot[lengthTot > min(range) & lengthTot < max(range)]
  tabl <- table(lengthTot)
  if(length(tabl) > 8) {
    tabl <- sort(tabl, decreasing = TRUE)[1:8]
  }
  # tabl <- as.table(tabl[tabl > cutoff * sum(tabl)])
  names(attributes(tabl)$dimnames) <- "AlleleLength"
  tabl <- tabl[order(as.numeric(names(tabl)))]
  return(tabl)
}

MostLikely <- function(results) {
  AlleleCandidate1 <- results[which.max(names(results))]
  resultsShort <- results[results != AlleleCandidate1]
  AlleleCandidate2 <- resultsShort[which.max(names(resultsShort))]
  results_sorted <- sort(results)
  max_peak <- results_sorted[length(results_sorted) - 0:1]
  AlleleCandidate3 <- max_peak[1]
  AlleleCandidate4 <- max_peak[2]
  AC <- list(AC1 = AlleleCandidate1[1], AC2 = AlleleCandidate2, 
             AC3 = AlleleCandidate3, AC4 = AlleleCandidate4)
  GT <- c(AC$AC1, AC$AC3)
  names(rev(GT))
}

PlotSample <- function(sampleName, sampleTab, sampleML) {
  barplot(sampleTab,
          col = BarCol(sampleTab, highlight = sampleML),
          main = sampleName,
          xlab = paste(sampleML, collapse = "/"))
}

PlotAllSamples <- function(loadedFiles, writeTo, primerInfo) {
  pdf(writeTo)
  outputDF <- data.frame()
  for(i in loadedFiles) {
    if(length(i[[1]]) == 0) {
      barplot(c(0, 0),
              main = unlist(i[2]),
              xlab = "0/0")
    } else {
      primerName <- unlist(strsplit(unlist(i[2]), "_"))[1]
      primerFR <- primerInfo[, colnames(primerInfo) == primerName]
      tabSample <- FragLength2(sample = i[[1]], primerF = primerFR[1], primerR = primerFR[2], range = 80:250)
      mlSample <- genotypefromNGS::MostLikely(results = tabSample)
      PlotSample(sampleName = unlist(i[2]), sampleTab = tabSample, sampleML = mlSample)
      outputDF <- rbind(outputDF, c(unlist(i[2]), mlSample))
    }
  }
  dev.off()
  return(outputDF)
}

test1 <- PlotAllSamples(loadedFiles = loadedFastqFiles, writeTo = "./Plots/results_reindeer_gt.pdf", primerInfo = primerDF)
test2 <- t(test1)

# sortera bort inkorrekta sekvenser av samma längd, längre distans än sekvensfel
