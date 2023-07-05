library(ShortRead)

primerDF <- data.frame(NVHRT03 = c("TGGAGAGCTGAGTATGAAAG", "CTTTTAGGTAGCTGCATTTCT"), 
                       NVHRT16 = c("ATTCTAAGCCCAAATAATCTT", "AAGACACAGACCCCTTAGA"), 
                       NVHRT73 = c("CTTGCCCATTTAGTGTTTTCT", "CTCCTATTCAATGACACGCA"),
                       NVHRT31 = c("CATCCCAAAGTTTACAGCAG", "GTGTCACTCAGGCCATAAAAA"),
                       BM6506 = c("GCACCTGGTAAAGAGATGGC", "GTGCCATGCTCAAGTTGCT"),
                       RT9 = c("TGAAGTTTAATTTCCACTCT", "ATGTGGGATGAAAGTGACTG"),
                       NVHRT01 = c("GCAGTCTTCCCCTTTCTT", "TAGTGTCCAACTCTGCAATC"),
                       OheQ = c("AGACCTGATTACAATGTGTCAGTGAAGGTCTTC", "CTAGATGGTTGCCTGGATGGGTCCATC"),
                       DeerC89 = c("AGAGCCTCGTCTTTTCATTC", "TTAGACAAGCAAGCAGCAAA"),
                       NVHRT66 = c("GCAGAGTCCGTGGGATTG", "ATAAGCCAAGCTGCCTCCAA"),
                       BM4513 = c("GCGCAAGTTTCCTCATGC", "GGGTGATGTACTGAATTGCTGA"),
                       NVHRT48 = c("CGTGAATCTTAACCAGGTCT", "GTTTCTAAATGAAGCTGACC"),
                       RT27 = c("CCAAAGACCCAACAGATG", "AATGCTTTTGCTGTGTTACAA"),
                       RT30 = c("CACTTGGCTTTTGGACTTA", "AGTGTGCATACATACACCAG"),
                       RT7 = c("CCTGTTCTACTCTTCTTCTC", "AACCAGTGCCCGTGAAAAGT"),
                       RT1 = c("TGCCTTCTTTCATCCAACAA", "GTAAAGAGGATGGGAAGATG"),
                       OarFCB193 = c("TTCATCTCAGACTGGGATTCAGAAAGGC", "GGGATGCAGGAGGGTTATTTCCAAGC"),
                       MAF46 = c("CACCATGGCCACCTGGAATCAGG", "GTGGTACTGTGCCTTATAGGGTATTT"))

inputFiles <- list.files("./filtered2/", pattern = "*.fastq.gz")

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
  loadedFile <- readFastq(dirPath = "./filtered2/", pattern = namePattern)
  print(paste("Searching for", namePattern))
  print(list.files("./filtered2/", pattern = namePattern))
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

test1 <- PlotAllSamples(loadedFiles = loadedFastqFiles, writeTo = "./plots/results_reindeer_gt.pdf", primerInfo = primerDF)
test2 <- t(test1)

# sortera bort inkorrekta sekvenser av samma längd, längre distans än sekvensfel