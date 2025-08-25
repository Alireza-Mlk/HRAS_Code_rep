################################
#
# R code used for GSEA using the fgsea package. Uses .tsv event files from rmats_pipeline.py
#
# RI_Events can be replaced by the event file to be analysed
#
# Requirements: dplyr, msigdbr, fgsea, tidyverse
#
library(dplyr)
library(msigdbr)
library(fgsea)
library(tidyverse)


frames_lis <- RI_Events



# Remove events from unknown genes
frames_list <- RI_Events_no_na %>%
  mutate(geneSymbol = (case_when(is.na(geneSymbol) ~ GeneID, !is.na(geneSymbol) ~ geneSymbol)))

# Function to calculate sample counts
get_counts <- function (x) {
  x <- strsplit(x, "-")
  x <- lapply(x, as.numeric)
  xsum <- lapply(x, sum)
  return(unlist(xsum))
  
}

significant_events_fg <- frames_list %>%

# Filter by FDR, Inclusion level difference, and sample counts. 
  mutate(sample1_count = get_counts(IJC_SAMPLE_1) + get_counts(SJC_SAMPLE_1),
         sample2_count = get_counts(IJC_SAMPLE_2) + get_counts(SJC_SAMPLE_2)) %>%
  filter( sample1_count > 50 & sample2_count > 50) %>%
  filter(TA == "T0" & TB != "T0VC")

# Replace 2.2e-16 p-value with zero (smallest number from rMATs)
  significant_events_fg$PValue[significant_events_fg$PValue == 0] <-  2.2e-16

# Separate events by second time points
  T1_fg <- filter(significant_events_fg, TB == "T1")
  T3_fg <- filter(significant_events_fg, TB == "T3")
  T5_fg <- filter(significant_events_fg, TB == "T5")
  T7_fg <- filter(significant_events_fg, TB == "T7")
  T1M_fg <- filter(significant_events_fg, TB == "T1M")
  T3M_fg <- filter(significant_events_fg, TB == "T3M")
  T5M_fg <- filter(significant_events_fg, TB == "T5M")
  T7M_fg <- filter(significant_events_fg, TB == "T7M")

# For each time point, rank events per gene by p-value and Inclusion level difference
# Remove duplicates and genes with only 0 inclusion level differences
  ### T1 ###
    rank_file <- T1_fg %>%

    mutate(T_z=paste(TA,"-",TB)) %>% 
    group_by(geneSymbol) %>% 
    filter(PValue == min(PValue)) %>% 
    filter(abs(IncLevelDifference) == max(abs(IncLevelDifference))) %>% 
    mutate(rank=IncLevelDifference) %>%
    ungroup() %>%
    select(geneSymbol, rank)
   
   rank_file <- rank_file[!duplicated(rank_file), ]
   rank_file <- group_by(rank_file, geneSymbol) %>%
     filter(mean(rank) != 0) %>%
     ungroup()

   rownms <- rank_file$geneSymbol
   rank_file <- select(rank_file, rank)
   ranks_1 <- rank_file$rank
   names(ranks_1) <- rownms
   
   
  ### T3 ###
   
   
   rank_file <- T3_fg %>%
     
     mutate(T_z=paste(TA,"-",TB)) %>% 
     group_by(geneSymbol) %>% 
     filter(PValue == min(PValue)) %>% 
     filter(abs(IncLevelDifference) == max(abs(IncLevelDifference))) %>% 
     mutate(rank=IncLevelDifference) %>%
     ungroup() %>%
     select(geneSymbol, rank)
   
   rank_file <- rank_file[!duplicated(rank_file), ] 
   rank_file <- group_by(rank_file, geneSymbol) %>%
    filter(mean(rank) != 0) %>%
     ungroup()
   
   rownms <- rank_file$geneSymbol
   rank_file <- select(rank_file, rank)
   ranks_3 <- rank_file$rank
   names(ranks_3) <- rownms
  
  ### T5 ###
   rank_file <- T5_fg %>%
     
     mutate(T_z=paste(TA,"-",TB)) %>% 
     group_by(geneSymbol) %>% 
     filter(PValue == min(PValue)) %>% 
     filter(abs(IncLevelDifference) == max(abs(IncLevelDifference))) %>% 
     mutate(rank=IncLevelDifference) %>%
     ungroup() %>%
     select(geneSymbol, rank)
   
   rank_file <- rank_file[!duplicated(rank_file), ]
   rank_file <- group_by(rank_file, geneSymbol) %>%
   filter(mean(rank) != 0) %>%
     ungroup()
   
   rownms <- rank_file$geneSymbol
   rank_file <- select(rank_file, rank)
   ranks_5 <- rank_file$rank
   names(ranks_5) <- rownms
   
   ### T7 ###
   rank_file <- T7_fg %>%
     
     mutate(T_z=paste(TA,"-",TB)) %>% 
     group_by(geneSymbol) %>% 
     filter(PValue == min(PValue)) %>% 
     filter(abs(IncLevelDifference) == max(abs(IncLevelDifference))) %>% 
     mutate(rank=IncLevelDifference) %>%
     ungroup() %>%
     select(geneSymbol, rank)
   
   rank_file <- rank_file[!duplicated(rank_file), ]
   rank_file <- group_by(rank_file, geneSymbol) %>%
   filter(mean(rank) != 0) %>%
     ungroup()
   
   rownms <- rank_file$geneSymbol
   rank_file <- select(rank_file, rank)
   ranks_7 <- rank_file$rank
   names(ranks_7) <- rownms
   
   ### T1M ###
   
   rank_file <- T1M_fg %>%
     
     mutate(T_z=paste(TA,"-",TB)) %>% 
     group_by(geneSymbol) %>% 
     filter(PValue == min(PValue)) %>% 
     filter(abs(IncLevelDifference) == max(abs(IncLevelDifference))) %>% 
     mutate(rank=IncLevelDifference) %>%
     ungroup() %>%
     select(geneSymbol, rank)
   
   rank_file <- rank_file[!duplicated(rank_file), ] 
   rank_file <- group_by(rank_file, geneSymbol) %>%
   filter(mean(rank) != 0) %>%
     ungroup()
   
   rownms <- rank_file$geneSymbol
   rank_file <- select(rank_file, rank)
   ranks_1M <- rank_file$rank
   names(ranks_1M) <- rownms
   
   ### T3M ###
   
   rank_file <- T3M_fg %>%
     
     mutate(T_z=paste(TA,"-",TB)) %>% 
     group_by(geneSymbol) %>% 
     filter(PValue == min(PValue)) %>% 
     filter(abs(IncLevelDifference) == max(abs(IncLevelDifference))) %>% 
     mutate(rank=IncLevelDifference) %>%
     ungroup() %>%
     select(geneSymbol, rank)
   
   rank_file <- rank_file[!duplicated(rank_file), ] 
   rank_file <- group_by(rank_file, geneSymbol) %>%
   filter(mean(rank) != 0) %>%
     ungroup()
   
   rownms <- rank_file$geneSymbol
   rank_file <- select(rank_file, rank)
   ranks_3M <- rank_file$rank
   names(ranks_3M) <- rownms
   
   ### T5M ###
   
   rank_file <- T5M_fg %>%
     
     mutate(T_z=paste(TA,"-",TB)) %>% 
     group_by(geneSymbol) %>% 
     filter(PValue == min(PValue)) %>% 
     filter(abs(IncLevelDifference) == max(abs(IncLevelDifference))) %>% 
     mutate(rank=IncLevelDifference) %>%
     ungroup() %>%
     select(geneSymbol, rank)
   
   rank_file <- rank_file[!duplicated(rank_file), ] 
   rank_file <- group_by(rank_file, geneSymbol) %>%
   filter(mean(rank) != 0) %>%
     ungroup()
   
   rownms <- rank_file$geneSymbol
   rank_file <- select(rank_file, rank)
   ranks_5M <- rank_file$rank
   names(ranks_5M) <- rownms
   
   ### T7M ###
   
   rank_file <- T7M_fg %>%
     
     mutate(T_z=paste(TA,"-",TB)) %>% 
     group_by(geneSymbol) %>% 
     filter(PValue == min(PValue)) %>% 
     filter(abs(IncLevelDifference) == max(abs(IncLevelDifference))) %>% 
     mutate(rank=IncLevelDifference) %>%
     ungroup() %>%
     select(geneSymbol, rank)
   
   rank_file <- rank_file[!duplicated(rank_file), ] 
   rank_file <- group_by(rank_file, geneSymbol) %>%
   filter(mean(rank) != 0) %>%
     ungroup()
   
   rownms <- rank_file$geneSymbol
   rank_file <- select(rank_file, rank)
   ranks_7M <- rank_file$rank
   names(ranks_7M) <- rownms

# Select geneset from Msigdb (C4: Computational, C5: Ontology, C6: Oncogenic)
   gmt_file <- msigdbr(species = "human", collection = "C6")
   gmt.file <- split(as.character(gmt_file$gene_symbol), gmt_file$gs_name)

   
# Run GSEA per timepoint 
  SE_T1_C6_GSEA <- fgsea(pathways = gmt.file, 
                     stats    = ranks_1,
                     eps      = 0,
                     minSize  = 10,
                     maxSize  = 100)
  
   SE_T3_C6_GSEA <- fgsea(pathways = gmt.file, 
                              stats    = ranks_3,
                              eps      = 0,
                              minSize  = 10,
                              maxSize  = 100)

   SE_T5_C6_GSEA <- fgsea(pathways = gmt.file, 
                              stats    = ranks_5,
                              eps      = 0,
                              minSize  = 10,
                              maxSize  = 100)
   
   SE_T7_C6_GSEA <- fgsea(pathways = gmt.file, 
                              stats    = ranks_7,
                              eps      = 0,
                              minSize  = 10,
                              maxSize  = 100)
   
   SE_T1M_C6_GSEA <- fgsea(pathways = gmt.file, 
                              stats    = ranks_1M,
                              eps      = 0,
                              minSize  = 10,
                              maxSize  = 100)
   
   SE_T3M_C6_GSEA <- fgsea(pathways = gmt.file, 
                               stats    = ranks_3M,
                               eps      = 0,
                               minSize  = 10,
                               maxSize  = 100)
   
   SE_T5M_C6_GSEA <- fgsea(pathways = gmt.file, 
                               stats    = ranks_5M,
                               eps      = 0,
                               minSize  = 10,
                               maxSize  = 100)

   SE_T7M_C6_GSEA <- fgsea(pathways = gmt.file, 
                               stats    = ranks_7M,
                               eps      = 0,
                               minSize  = 10,
                               maxSize  = 100)   
