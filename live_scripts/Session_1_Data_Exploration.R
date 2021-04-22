library(tximport)
library(DESeq2)
library(tidyverse)

# sample meta data

sampleinfo <- read_tsv("data/samplesheet.tsv")

# read in the count

files <- str_c("salmon/", sampleinfo$SampleName, "/quant.sf")
files
files <- set_names(files, sampleinfo$SampleName)
files

tx2gene <- read_tsv("references/tx2gene.tsv")
tx2gene

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
str(txi)

saveRDS(txi, file = "salmon_outputs/txi.rds")

# Create a raw counts matrix

rawCounts <- round(txi$counts, 0)
head(rawCounts)
dim(rawCounts)

keep <- rowSums(rawCounts) > 5
table(keep)
head(keep)

filtCounts <- rawCounts[keep, ]
dim(filtCounts)

# Distribution of the counts

summary(filtCounts)

boxplot(filtCounts, main = "Raw Counts", las=2)

plot(rowMeans(filtCounts), rowSds(filtCounts), xlim=c(0, 10000), ylim=c(0,5000))

# Data transformation

logcounts <- log2(filtCounts + 1)

boxplot(logcounts, main = "Log Counts", las=2)

plot(rowMeans(logcounts), rowSds(logcounts))

# VST - variance stabilising transformation

vstCounts <- vst(filtCounts)

boxplot(vstCounts, main = "VST Counts", las=2)

plot(rowMeans(vstCounts), rowSds(vstCounts))

# rlog - Exercise 2

?rlog

rlogcounts <- rlog(filtCounts)
head(rlogcounts)

boxplot(rlogcounts, main = "regularised log Counts", las=2)

## Principal Component Analysis

library(ggfortify)

pcDat <- prcomp(t(rlogcounts))
str(pcDat)

autoplot(pcDat)

### add some colour and shapes

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5)

# Exercise 3

?autoplot.prcomp

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5,
         x = 2,
         y = 3)

## Discussion on the PCA

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5)

### Which are the mislabelled samples

library(ggrepel)

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5) +
  geom_text_repel(aes(x = PC1, y = PC2, label=SampleName),
                  box.padding = 0.8)

sampleinfo <- mutate(sampleinfo, Status = case_when(
  SampleName == "SRR7657882" ~ "Uninfected",
  SampleName == "SRR7657873" ~ "Infected",
  TRUE ~ Status
))

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5) 

write_tsv(sampleinfo, "results/SampleInfo_corrected.tsv")


