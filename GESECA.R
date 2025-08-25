library(dplyr)
library(msigdbr)
library(fgsea)
library(tidyverse)
library(preprocessCore)
library(ggplot2)

na_data <- SE_Events
filter_data <- filter(na_data, is.na(geneSymbol))
nfilter_data <- filter(na_data, !is.na(geneSymbol))
ensmbl_G <- ensmbl_gtf
for (i in 1:nrow(filter_data)) {
  # if (is.na(filter_data$geneSymbol[i])) {
  chr = filter_data$chr[i]
  # print(i)
  # print("Loop 1")
  slct_gtf <- filter(ensmbl_G, seqid == chr)
  for (j in 1:nrow(slct_gtf)) {
    # print(j)
    # print("Loop 2")
    if (filter_data$exonStart_0base[i] >= slct_gtf$start[j] && filter_data$exonEnd[i] <= slct_gtf$end[j]) {
      filter_data$geneSymbol[i] <- slct_gtf$gene_name[j]
      # print("If Statement here!!!!")
      break 
    }
  }
}
frames_list <- rbind(filter_data, nfilter_data)



frames_list <- RI_Events_no_na %>%
  mutate(geneSymbol = (case_when(is.na(geneSymbol) ~ GeneID, !is.na(geneSymbol) ~ geneSymbol)))

get_counts <- function (x) {
  x <- strsplit(x, "-")
  x <- lapply(x, as.numeric)
  xsum <- lapply(x, sum)
  return(unlist(xsum))
  
}

significant_events_fgn <- mutate(frames_list, T_z=paste(TA,"-",TB, sep = "")) 
  significant_events_fgn <- significant_events_fgn %>%
    filter(!grepl("M", TB) & TA == "T0" & TB != "T0VC") %>%
    # filter(TA == "T0" & TB != "T0VC" ) %>%
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
  significant_events_fgn <- significant_events_fgn[!duplicated(significant_events_fgn), ]
  
  significant_events_fgn <- filter(significant_events_fgn, mean(IncLevelDifference) != 0) %>%
    ungroup()

  fgn_m <- pivot_wider(significant_events_fgn, id_cols = geneSymbol, names_from = T_z, values_from = IncLevelDifference)
  fgn_m$`T0-T1`[is.na(fgn_m$`T0-T1`)] <- 0
  fgn_m$`T0-T3`[is.na(fgn_m$`T0-T3`)] <- 0
  fgn_m$`T0-T5`[is.na(fgn_m$`T0-T5`)] <- 0
  fgn_m$`T0-T7`[is.na(fgn_m$`T0-T7`)] <- 0
  
  
  rank_m <- select(fgn_m, -geneSymbol) 
  rownames(rank_m) = fgn_m$geneSymbol
  rank_m <- as.matrix(rank_m)

  rank_mq <-  normalize.quantiles(rank_m)
  rownames(rank_mq) <- rownames(rank_m)
  colnames(rank_mq) <- colnames(rank_m)
  
  rank_mq <- as.matrix(rank_mq)
  
  gmt_m<- msigdbr(species = "human", collection = "C4")
  gmt.m <- split(as.character(gmt_m$gene_symbol), gmt_m$gs_name)
  

  geseca_C4 <- geseca(gmt.m, rank_mq, minSize = 15, maxSize = 500)
  ### VISUALISATION-2
  plotGesecaTable(geseca_C4[order(geseca_C4$padj, decreasing = FALSE)] |> head(15), gmt.m, E=rank_mq)
  
  
  
  
    ### VISUALISATION-1
  

  
  plotCoregulationProfile(gmt.m[["GSE15330_MEGAKARYOCYTE_ERYTHROID_VS_GRANULOCYTE_MONOCYTE_PROGENITOR_DN"]], rank_mq, conditions = colnames(rank_mq))
  

  
  
  
  
  
  
   # fgn_file <- significant_events_fgn %>%

    # mutate(T_z=paste(TA,"-",TB)) %>%
    # group_by(geneSymbol) %>%
    # filter(PValue <= abs(maximum(PValue)))
    # mutate(rank=IncLevelDifference*-log10(PValue)) %>%
    # select(geneSymbol, rank)

    # group_by(geneSymbol) %>% 
    # filter(PValue == min(PValue)) %>% 
    # filter(abs(IncLevelDifference) == max(abs(IncLevelDifference))) %>% 
    # ungroup() %>%
    # select(geneSymbol, IncLevelDifference, T_z)
  # fgn_file <- fgn_file[!duplicated(fgn_file), ]
  