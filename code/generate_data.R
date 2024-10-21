library(tidyverse)
library(tidylog)
library(readxl)
library(qiime2R)
library(ggprism)

theme_set(theme_prism())
unknown_samples <- c('AB001', 'AE1647', 'AE1364', 'AE1402')

# Load in read counts
read_counts <- read_tsv('data/240711_merged/data/sequences/per-sample-fastq-counts_trimmed.tsv') %>%
  rename(trimmed_forward = `forward sequence count`,
         trimmed_reverse = `reverse sequence count`,
         sampleid = `sample ID`) %>%
  filter(!sampleid %in% unknown_samples)

# # Dropping samples with < 10,000 reads
samples_to_drop <- read_counts %>%
  filter(trimmed_forward < 1e4)


table <- read_q2biom('data/240711_merged/data/merged/feature-table_taxa.biom') %>%
  as_tibble(rownames = "featureid") %>%
  select(-any_of(c(unknown_samples, sample_to_drop$sampleid)))

breed_list <- read_csv('data/breed_list.csv')

acute_gi_list <- read_excel('data/240711_merged/data/acute_gi_classifications.xlsx') %>%
  select(-diagnosis)

metadata <- read_xlsx('data/240711_merged/metadata_merged.xlsx', col_types = 'text') %>%
  janitor::clean_names() %>%
  separate(sampleid, into = "sampleid", sep = "_") %>% # drop plate location
  select(-c(source_material_id, plate_id, sra_number)) %>% # drop uninformative metadata
  mutate(
    dob = as.Date(as.double(dob), origin = "1899-12-30"),
    collection_date = as.Date(as.double(collection_date), origin = "1899-12-30"),
    age_of_collection = difftime(collection_date, dob, units = 'days'),
    month_norm = month(collection_date) + (day(collection_date) / days_in_month(collection_date)),
    sample_order = as.double(sample_order),
    visit = as.double(visit),
    disease_category = str_to_title(disease_category),
    feed_hay = as.double(feed_hay),
    feed_pasture = as.double(feed_pasture),
    feed_grain = as.double(feed_grain),
    feed_nursing = as.double(feed_nursing),
  ) %>%
  separate(geo_loc_name, into = c('country', 'university'), sep = ': ') %>%
  group_by(horse, visit) %>%
  arrange(collection_date) %>%
  mutate(sample_age = as.double(1 + age_of_collection - first(age_of_collection))) %>%
  ungroup() %>% 
  # filter(!sampleid %in%  samples_to_drop$sampleid) %>%
  # filter(!sampleid %in% unknown_samples) %>%
  # filter(sampleid %in% colnames(table)) %>%
  left_join(., breed_list, by = 'breed') %>%
  left_join(., acute_gi_list, by = c('horse', 'visit'))

# setdiff(metadata$sampleid, colnames(table)) # all metadata sampleids in table
# setdiff(colnames(table), metadata$sampleid) # all table columns in metadata

metadata %>%
  write_rds('data/240711_merged/data/metadata.Rds')
# table %>%
#   write_rds('data/240711_merged/data/dada2_table.Rds')

# taxonomy <-  read_tsv('data/240711_merged/data/taxonomy/taxonomy.tsv') %>%
#     rename(featureid = `Feature ID`) %>%
#   select(featureid, Taxon) %>%
#   separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
#            sep = '; ') %>%
#   pivot_longer(-featureid, names_to = 'level',
#                values_to = 'taxon') %>%
#   mutate(taxon = substr(taxon, 4, nchar(taxon))) %>%
#   group_by(featureid) %>%
#   fill(taxon) %>%
#   ungroup()
# taxonomy %>%
#   write_rds('data/240711_merged/data/taxonomy.Rds')

metadata <- read_rds('data/240711_merged/data/metadata.Rds')
table <- read_rds('data/240711_merged/data/dada2_table.Rds')
taxonomy <- read_rds('data/240711_merged/data/taxonomy.Rds')
