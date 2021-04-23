# load libraries
library(DESeq2)
library(tidyverse)

# load data
txi <- readRDS( 'RObjects/txi.rds' )
sampleinfo <- read_tsv( 'data/samplesheet_corrected.tsv', col_types = c('cccc') )

class(txi)
names(txi)
head(txi$counts)
head(txi$length)
head(sampleinfo)

# sample order
all(colnames(txi$counts) == sampleinfo$SampleName )

# Design
simple.model <- as.formula( ~Status)
model.matrix( simple.model, data=sampleinfo)
as.factor(sampleinfo$Status)

# relevel
sampleinfo <- mutate( sampleinfo, Status = fct_relevel(Status, c('Uninfected', 'Infected')))
model.matrix( simple.model, data=sampleinfo)

# build DDS object
ddsObj.raw <- DESeqDataSetFromTximport( txi = txi,
                                        colData = sampleinfo,
                                        design = simple.model)

# filter
keep <- rowSums( counts(ddsObj.raw)) > 5
dim(counts(ddsObj.raw))
ddsObj.filt <- ddsObj.raw[keep, ]

dim(counts(ddsObj.filt))

# estimate size factors
ddsObj <- estimateSizeFactors( ddsObj.filt )

normalizationFactors( ddsObj.filt)
head(normalizationFactors(ddsObj))

apply( normalizationFactors( ddsObj), 2, median)

# MA plot: raw data
logcounts <- log2( counts(ddsObj, normalized = FALSE) + 1)
head(logcounts)

limma::plotMA( logcounts, array = 5, ylim= c(-5,5))
abline(h=0, col='red')

# MA plot : Normalised dada
logcounts <- log2( counts(ddsObj, normalized = TRUE) + 1)
head(logcounts)

limma::plotMA( logcounts, array = 5, ylim= c(-5,5))
abline(h=0, col='blue')

# estimate disperssions
ddsObj <- estimateDispersions(ddsObj)

plotDispEsts(ddsObj)

# Var and dispr
# var = mean + disp + mean^2

# fit and test
ddsObj <- nbinomWaldTest(ddsObj)

# DEseq = estimateSizeFactors + estimateDispersions + nbinomWaldTest
ddsObj <- DESeq( ddsObj.filt)

# results table
results.simple <- results( ddsObj, alpha = 0.05)

# how many genes sig. DE
table( results.simple$padj < 0.05, useNA='always')

# p-val
table( results.simple$pvalue < 0.05, useNA='always')

# how manu genes are up regulated
sum(results.simple$padj < 0.05 & results.simple$log2FoldChange > 0, na.rm=T)
# down regulated
sum(results.simple$padj < 0.05 & results.simple$log2FoldChange < 0, na.rm=T)
########################################################################
# Staus = 4 exp
# TimePoint = 3 exp
# Status + TimePoint = 7

additive.model <- as.formula( ~TimePoint + Status )
model.matrix( additive.model, data=sampleinfo)

# power and parameter estimation
# var = SS / df
# df = n - no. parameters estimate

# build DDS object
ddsObj.raw <- DESeqDataSetFromTximport( txi = txi,
                                        colData = sampleinfo,
                                        design = additive.model)
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]
dim(counts(ddsObj.filt))

ddsObj <- DESeq( ddsObj.filt)
# results

results.additive <- results(ddsObj, alpha = 0.05)

sum(results.additive$padj < 0.05, na.rm=T)
sum(results.simple$padj < 0.05, na.rm=T)

resultsNames(ddsObj)

results.InfectedvUninfected <- results( ddsObj, name='Status_Infected_vs_Uninfected')

# topGenes
topGenesIvU <- as.data.frame(results.InfectedvUninfected) %>% 
  rownames_to_column( 'GeneID' ) %>% 
  top_n(100, wt=-padj)
  

dim(topGenesIvU)

# PCA
vastcounts <- vst( ddsObj.raw, blind = TRUE)
design(ddsObj.raw)
plotPCA(vastcounts, intgroup = c( 'Status', 'TimePoint'))

# compare models
ddsObj.LRT <- DESeq( ddsObj, test='LRT', reduced = simple.model)

sum( results.Additive_v_Simple$padj < 0.05, na.rm=TRUE)

results.Additive_v_Simple <- results( ddsObj.LRT)
results.Additive_v_Simple

# Interaction model

interaction.model <- as.formula( ~TimePoint + Status + TimePoint:Status)
interaction.model <- as.formula( ~TimePoint * Status )
model.matrix(interaction.model, data=sampleinfo)


# build DDS object
ddsObj.raw <- DESeqDataSetFromTximport( txi = txi,
                                        colData = sampleinfo,
                                        design = interaction.model)
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]
dim(counts(ddsObj.filt))

ddsObj.interaction <- DESeq( ddsObj.filt)

# Interaction Vs additive 
ddsObj.LRT <- DESeq( ddsObj.interaction, test='LRT', reduced = additive.model)

# extract the results
results.Interaction_v_Additive <- results( ddsObj.LRT)

table(results.Interaction_v_Additive$padj < 0.05, useNA='always')

resultsNames(ddsObj.interaction)

results.interaction.11 <- results( ddsObj.interaction,
                                   name = 'Status_Infected_vs_Uninfected',
                                   alpha = 0.05
                                   )

sum(results.interaction.11$padj < 0.05, na.rm=T)

# Inf Vs Unin 33 days post infec
results.interaction.33 <- results(ddsObj.interaction,
                                  contrast = list( c('Status_Infected_vs_Uninfected', 
                                                     'TimePointd33.StatusInfected')
                                              ),
                                  alpha = 0.05
                                  )

sum(results.interaction.33$padj < 0.05, na.rm=T)

# d33 vs d11 in Uninfectd
results.d33_v_d11_Uninfected <- results(ddsObj.interaction,
                                        name='TimePoint_d33_vs_d11',
                                        alpha = 0.05
                                        )

sum(results.d33_v_d11_Uninfected$padj < 0.05, na.rm=T)

# d33 vs d11 in Infected mice
results.interaction.33 <- results(ddsObj.interaction,
                                  contrast = list( c('TimePoint_d33_vs_d11', 
                                                     'TimePointd33.StatusInfected')
                                  ),
                                  alpha = 0.05
)

sum(results.interaction.33$padj < 0.05, na.rm=T)

# save robjects
saveRDS( ddsObj.interaction, 'results/DEseqDataSet.interaction.rds')

saveRDS(results.d33_v_d11_Uninfected, 'results/DEseqResults.interaction_d11.rds')
saveRDS(results.interaction.33, 'results/DEseqResults.interaction_d33.rds')


