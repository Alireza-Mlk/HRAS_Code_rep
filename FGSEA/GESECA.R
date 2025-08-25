################################
#
# R code used for GESECA using the fgsea package. Uses .tsv event files from rmats_pipeline.py
#
# SE_Events can be replaced by the event file to be analysed
#
# Requirements: dplyr, msigdbr, fgsea, tidyverse, preprocessCore, ggplot2
#
library(dplyr)
library(msigdbr)
library(fgsea)
library(tidyverse)
library(preprocessCore)
library(ggplot2)

# Removes events from unknown genes
frames_list <- SE_Events %>%
  mutate(geneSymbol = (case_when(is.na(geneSymbol) ~ GeneID, !is.na(geneSymbol) ~ geneSymbol)))

# Function to calculate sample counts
get_counts <- function (x) {
  x <- strsplit(x, "-")
  x <- lapply(x, as.numeric)
  xsum <- lapply(x, sum)
  return(unlist(xsum))
  
}

# Filter specific timepoints
significant_events_fgn <- mutate(frames_list, T_z=paste(TA,"-",TB, sep = "")) 
  significant_events_fgn <- significant_events_fgn %>%
    filter(!grepl("M", TB) & TA == "T0" & TB != "T0VC") %>%
# Filter by FDR, Inclusion level difference, and sample counts. And filter the most significant event for each gene per timepoint.
  group_by(geneSymbol) %>%
  filter(max(abs(IncLevelDifference)) > 0.1 & min(FDR) < 0.05) %>%
  mutate(sample1_count = get_counts(IJC_SAMPLE_1) + get_counts(SJC_SAMPLE_1),
         sample2_count = get_counts(IJC_SAMPLE_2) + get_counts(SJC_SAMPLE_2)) %>%
  filter( sample1_count > 50 & sample2_count > 50) %>%
  ungroup() %>%
  group_by(T_z) %>%
  group_by(geneSymbol, .add = TRUE) %>%
  filter(PValue == min(PValue)) %>% 
  filter(abs(IncLevelDifference) == max(abs(IncLevelDifference))) %>%
  select(geneSymbol, IncLevelDifference, T_z)

# Remove duplicates
  significant_events_fgn <- significant_events_fgn[!duplicated(significant_events_fgn), ]

# Inclusion level difference is set to 0 for mirrored events with equally - and + inclusion level differences for the same gene  

  significant_events_fgn <- filter(significant_events_fgn, mean(IncLevelDifference) != 0) %>%
    ungroup()

# Replace NA with 0
  fgn_m <- pivot_wider(significant_events_fgn, id_cols = geneSymbol, names_from = T_z, values_from = IncLevelDifference)
  fgn_m$`T0-T1`[is.na(fgn_m$`T0-T1`)] <- 0
  fgn_m$`T0-T3`[is.na(fgn_m$`T0-T3`)] <- 0
  fgn_m$`T0-T5`[is.na(fgn_m$`T0-T5`)] <- 0
  fgn_m$`T0-T7`[is.na(fgn_m$`T0-T7`)] <- 0
  
# Normalise results and put into a matrix
  rank_m <- select(fgn_m, -geneSymbol) 
  rownames(rank_m) = fgn_m$geneSymbol
  rank_m <- as.matrix(rank_m)

  rank_mq <-  normalize.quantiles(rank_m)
  rownames(rank_mq) <- rownames(rank_m)
  colnames(rank_mq) <- colnames(rank_m)
  
  rank_mq <- as.matrix(rank_mq)

# Select geneset from Msigdb (C4: Computational, C5: Ontology, C6: Oncogenic)
  gmt_m<- msigdbr(species = "human", collection = "C4")
  gmt.m <- split(as.character(gmt_m$gene_symbol), gmt_m$gs_name)
  
# Rung GESECA
  geseca_C4 <- geseca(gmt.m, rank_mq, minSize = 15, maxSize = 500)
  ### VISUALISATION-1 TOP 15 pathways
  plotGesecaTable(geseca_C4[order(geseca_C4$padj, decreasing = FALSE)] |> head(15), gmt.m, E=rank_mq)
  
  
  
  
    ### VISUALISATION-2 specific pathway from visualisation-1
  
  plotCoregulationProfile(gmt.m[["GSE15330_MEGAKARYOCYTE_ERYTHROID_VS_GRANULOCYTE_MONOCYTE_PROGENITOR_DN"]], rank_mq, conditions = colnames(rank_mq))
