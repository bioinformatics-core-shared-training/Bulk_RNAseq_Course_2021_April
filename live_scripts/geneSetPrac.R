# over-representation analysis
#-----------------------------------

# toy example
# 20.000 genes
# 100 DEGs
# gene set A: 2000 genes
# out of these 2000 genes in gene set A 20 are also DE.

# contingency table

contingengyTable <- data.frame(
  NoDiffExp = c(17920, 1980),
  YesDiffExp = c(80, 20)
)
contingengyTable
row.names(contingengyTable) <- c("NoPathway", "YesPathway")
contingengyTable


# perform test
fisher.test(contingengyTable, alternative="greater")

# clusterProfiler to perform test on many gene sets at once
# KEGG

library(clusterProfiler)
library(tidyverse)

# load data: shrunken lfc

shrink.d11 <- readRDS("RObjects/Shrunk_Results.d11.rds")

# select DEGs

sigGenes <- shrink.d11 %>%
  drop_na(Entrez, FDR) %>%
  filter(FDR < 0.05 & abs(logFC) > 1) %>%
  pull(Entrez)
sigGenes %>% head()

kk <- enrichKEGG(gene=sigGenes, organism='mmu')
class(kk)
head(kk, n=10)

# visualise pathway

browseKEGG(kk, "mmu04612")

# add logFC data
# using pathview

library(pathview)
logFC <- shrink.d11$logFC
head(logFC)
names(logFC)  <- shrink.d11$Entrez
head(logFC)

pathview(gene.data = logFC,
         pathway.id = "mmu04612",
         species = "mmu",
         limit = list(gene=20, cpd=1))

# exercise 1: run pathview with 'mmu04659', with genes selected at FDR < 0.01

# select genes:
logFC <- shrink.d11 %>%
  drop_na(FDR, Entrez) %>%
  filter(FDR < 0.01) %>%
  dplyr::select(Entrez, logFC) %>%
  deframe()
head(logFC)  

pathview(gene.data = logFC,
         pathway.id = 'mmu04659',
         species = 'mmu',
         limit = list(gene=20, cpd=1))

# (exercise 2 - same but with GO - skipped - see github)

# gene set enrichment analysis
#-----------------------------------

# GSEA and MSigDB

library(msigdbr)

# rank genes
# by logFC

rankGenes <- shrink.d11 %>%
  drop_na(Entrez) %>%
  mutate(rank = logFC) %>%
  arrange(-rank) %>%
  pull(rank, Entrez)
head(rankGenes)  
tail(rankGenes)  
  
# load pathways
# hallmark
# need 'term-to-gene' info: list of gene in each pathway

m_H_t2g <- msigdbr(species = 'Mus musculus',
                   category = 'H') %>%
  dplyr::select(gs_name, entrez_gene, gene_symbol)
head(m_H_t2g)  

# perform analysis

gseaRes <- GSEA(rankGenes,
                TERM2GENE = m_H_t2g,
                pvalueCutoff = 1,
                minGSSize = 5,
                maxGSSize = 500)
head(gseaRes)
gseaRes %>% as_tibble

format.e1 <- function(x) {sprintf("%.1e", x)}

gseaRes %>%
  arrange(desc(abs(NES))) %>%
  top_n(10, -p.adjust) %>%
  dplyr::select(-core_enrichment) %>%
  dplyr::select(-Description) %>%
  data.frame() %>%
  mutate(ES=formatC(enrichmentScore, digits=3)) %>%
  mutate(NES=formatC(NES, digits=3)) %>%
  modify_at(
    c("pvalue", "p.adjust", "qvalues"),
    format.e1) %>%
  DT::datatable(rownames = FALSE, options = list(dom='t'))

# draw plot
# HALLMARK_INFLAMMATORY_RESPONSE

topx <- match("HALLMARK_INFLAMMATORY_RESPONSE",
              data.frame(gseaRes)$ID)
topx
head(gseaRes$ID)

# plot ranked genes
gseaplot(gseaRes, geneSetID = topx, by = 'preranked')

# plot ruuning score
gseaplot(gseaRes, geneSetID = topx, by = 'runningScore')

# Exercise 3

# rank genes using significance of the test for diff exp.
# logFC and p value

rankGenes.e1 <- shrink.d11 %>%
  drop_na(Entrez, pvalue, logFC) %>%
  mutate(rank = -log10(pvalue)*sign(logFC)) %>%
  arrange(-rank) %>%
  pull(rank, Entrez)
head(rankGenes.e1)

gseaRes.e1 <- GSEA(rankGenes.e1,
                TERM2GENE = m_H_t2g,
                pvalueCutoff = 1,
                minGSSize = 5,
                maxGSSize = 500)
head(gseaRes.e1)
gseaRes.e1 %>% as_tibble

gseaRes.e1 %>%
  arrange(desc(abs(NES))) %>%
  top_n(10, -p.adjust) %>%
  dplyr::select(-core_enrichment) %>%
  dplyr::select(-Description) %>%
  data.frame() %>%
  mutate(ES=formatC(enrichmentScore, digits=3)) %>%
  mutate(NES=formatC(NES, digits=3)) %>%
  modify_at(
    c("pvalue", "p.adjust", "qvalues"),
    format.e1) %>%
  DT::datatable(rownames = FALSE, options = list(dom='t'))
