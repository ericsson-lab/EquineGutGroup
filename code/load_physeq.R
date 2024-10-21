{
  source('code/load_data.R')
  library(tidyverse)
  library(phyloseq)
  library(microbiomeMarker)
  library(microViz)
  library(beepr)
}

tax_table <- taxonomy %>% 
  mutate(taxon = str_replace_all(taxon, '\\.', '_')) %>% 
  pivot_wider(names_from = 'level', values_from = 'taxon') %>% 
  filter(Kingdom == 'Bacteria') %>% 
  mutate(Family = paste(Phylum, Family, sep = '.')) %>%
  # rename(Domain = "Kingdom") %>% 
  column_to_rownames(var = 'featureid') %>% 
  as.matrix()

asv_table <- table %>% 
  column_to_rownames(var = 'featureid') %>% 
  as.matrix()

sample_metadata <- metadata %>% 
  group_by(horse, visit) %>% 
  mutate(samples = n(),
         SampleID = sampleid) %>%
  mutate(age = age_of_collection / 365) %>% 
  column_to_rownames(var='sampleid') %>% 
  as.data.frame()

colnames(sample_metadata)
ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(tax_table)
METADATA <- sample_data(sample_metadata) 

physeq <- phyloseq(ASV, TAX, METADATA)

physeq %>% 
  write_rds('data/240711_merged/physeq.Rds')


