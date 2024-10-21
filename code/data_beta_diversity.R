{
  library(vegan)
  library(ape)
  # source("code/load_data.R")
  
}

generate_dist <- function(table, distance = "bray"){
  
  table_t <- table %>% 
    column_to_rownames(var = "featureid") %>% 
    t()
  
  table.transformed <- table_t
  
  dist <- vegan::vegdist(table.transformed, method= distance,)
}

generate_pcoa <- function(dist, metadata_df) {
  pcoa_res <- ape::pcoa(dist, correction = "cailliez") 
  
  p_var <- (pcoa_res$values$Eigenvalues/pcoa_res$trace)*100
  
  pcoa_vectors <- pcoa_res$vectors %>% 
    as_tibble(rownames = "sampleid") %>% 
    select(sampleid, Axis.1, Axis.2, Axis.3, Axis.4, Axis.5)
  
  colnames(pcoa_vectors) <- c("sampleid", "PCo1", "PCo2", "PCo3", 
                              "PCo4", "PCo5")
  
  pcoa.metadata <- left_join(pcoa_vectors, metadata_df, by = "sampleid") 

  output <- list(pcoa.metadata, p_var)
  return(output)
}

# bc_dist <- generate_dist(table = table, distance = 'bray')
# bc_dist %>%
  # write_rds('data/bc_dist_raw.Rds')
# bc_dist %>% 
#   write_rds('data/bc_dist.Rds')

# beepr::beep()
