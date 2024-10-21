library(rstatix)
source('code/load_data.R')
source('code/data_beta_diversity.R')

theme_set(ggprism::theme_prism())

acute_gi_samples <- metadata %>% 
  filter(
    disease_category_acute_gi == 1 |
           disease_category_healthy == 1) %>% 
  filter(age_of_collection >= 4 * 365) %>% 
  filter(sample_age <= 7) %>% 
  mutate(acute_gi_subgroup = str_to_sentence(acute_gi_subgroup),
         disease_category = case_when(disease_category != 'Healthy' ~ 'Acute GI',
                                      TRUE ~ 'Healthy'),
         disease_category = factor(disease_category, levels = c('Healthy', 'Acute GI')),
         sample_age = as.double(sample_age)) %>% 
  filter(sampleid != 'JW1')
  

alpha_stats <- read_rds('data/240711_merged/data/alpha.stats.Rds') %>% 
  filter(sampleid %in% acute_gi_samples$sampleid) %>% 
  left_join(., acute_gi_samples)

### Metadata Summary
# Supplementary Figure 6a
acute_gi_samples %>% 
  select(disease_category, horse, sex) %>% 
  distinct() %>% 
  count(disease_category, sex) %>% 
  ggplot(aes(x = disease_category, y = n, fill = interaction(disease_category, sex))) +
  geom_col(stat = 'identity', color = 'black') +
  scale_fill_manual(values = c('blue', 'red', '#8080ff', '#ff8080')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  ylab('# of Patients') +
  theme(
    aspect.ratio = 2/1,
    axis.title.x = element_blank(),
    legend.position = 'none'
  )
ggsave('plots/acute_gi/sex_summary.png', width = 3, height = 4)

# Supplementary Figure 6b
acute_gi_samples %>% 
  select(disease_category, horse, country) %>% 
  distinct() %>% 
  ggplot(aes(x = fct_rev(fct_infreq(country)), fill = disease_category)) +
  geom_bar( position = position_dodge(), width = 0.5) +
  scale_fill_manual(values = c('blue', 'red'), name = 'Diagnosis') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = c(0, 150, 300),
                     limits = c(0, 300)) +
  ylab('# of Patients') +
  theme(
    aspect.ratio = 2/1,
    axis.title.y = element_blank(),
    legend.text = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold'),
    legend.position = c(0.65,0.15)
  ) +
  coord_flip()

ggsave('plots/acute_gi/country_summary.png', width = 4, height = 4)

# Supplementary Figure 6c
acute_gi_samples %>% 
  select(disease_category, horse, breed_type) %>% 
  distinct() %>% 
  ggplot(aes(x = fct_rev(fct_infreq(breed_type)), fill = disease_category)) +
  geom_bar( position = position_dodge(), width = 0.5) +
  scale_fill_manual(values = c('blue', 'red'), name = 'Diagnosis') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = c(0, 50, 100, 150),
                     limits = c(0, 150)) +
  ylab('# of Patients') +
  theme(
    aspect.ratio = 2/1,
    axis.title.y = element_blank(),
    legend.text = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold'),
    legend.position = c(0.65,0.15)
  ) +
  coord_flip()
ggsave('plots/acute_gi/breed_type_summary.png', width = 4, height = 4)

# Supplementary Figure 6d
acute_gi_samples %>% 
  ggplot(aes(x= age_of_collection/365, fill = fct_rev(disease_category))) +
  geom_histogram(position = 'identity', alpha = 0.5, 
                 data = acute_gi_samples %>% 
                   filter(disease_category == 'Acute GI')) +
  geom_histogram(position = 'identity', alpha = 0.7, 
                 data = acute_gi_samples %>% 
                   filter(disease_category == 'Healthy')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c('blue', 'red'), name = 'Diagnosis',
                    breaks = c('Healthy', 'Acute GI')) +
  xlab('Age at Collection (Years)') +
  ylab('# of Samples') +
  theme(
    aspect.ratio = 2/1,
    legend.text = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold'),
    legend.position = c(0.65,0.85)
  ) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 1))
  ) +
  coord_flip()
ggsave('plots/acute_gi/age_summary.png', width = 4, height = 4)



alpha_stats %>% 
  group_by(horse, visit) %>% 
  filter(n() >= 3) %>% 
  ggplot(aes(x = sample_age, y = chao1, 
             group = interaction(horse, visit), color = disease_category)) +
  geom_line(alpha = 0.3) +
  stat_summary(geom = 'line', fun = 'mean', linewidth = 1.5,
               aes(group = disease_category)) +
  scale_color_manual(values = c('blue', 'red'), name = 'Diagnosis') +
  scale_x_continuous(breaks = seq(1,7,1), 
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(0,3000), 
                     expand = expansion(mult = c(0, 0.05))) +
  xlab('Day') +
  ylab('Chao1 Index') +
  theme(
    aspect.ratio = 1/1.15,
    # legend.position = 'none'
  )
ggsave('plots/acute_gi/acute_gi_chao1_1wk.png', width = 5.25, height = 5)


alpha_stats %>% 
  group_by(horse, visit) %>% 
  filter(n() >= 3) %>% 
  ggplot(aes(x = sample_age, y = shannon, 
             group = interaction(horse, visit), color = disease_category)) +
  geom_line(alpha = 0.3) +
  stat_summary(geom = 'line', fun = 'mean', linewidth = 1.5,
               aes(group = disease_category)) +
  scale_color_manual(values = c('blue', 'red'), name = 'Diagnosis') +
  scale_x_continuous(breaks = seq(1,7,1), 
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(0,8), 
                     expand = expansion(mult = c(0, 0.05))) +
  xlab('Day') +
  ylab('Shannon Index') +
  theme(
    aspect.ratio = 1/1.5,
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16)
    # legend.position = 'none'
  )
ggsave('plots/acute_gi/acute_gi_shannon_1wk.png', width = 7, height = 5)
  


## Single Time point
alpha_stats %>%
  filter(sample_age == 1) %>%
  ggplot(aes(x = disease_category, y = chao1, color = disease_category)) +
  geom_point(position = position_jitter(0.2), size = 3, 
             alpha = 0.5) +
  stat_summary(geom = 'bar', fun = mean, fill = NA, 
               color = 'black', width = 0.4, linewidth = 1.5) +
  stat_summary(geom = 'errorbar', 
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               color = 'black', width = 0.2, linewidth = 1.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  # scale_x_discrete(labels = c('Healthy', 'Acute GI')) +
  scale_color_manual(values = c( 'blue', 'red')) +
  ylab('Chao1 Index') +
  theme(
    aspect.ratio = 2/1,
    legend.position = 'none',
    axis.title.x = element_blank()
  )

ggsave('plots/acute_gi/chao1_timepoint_1.png', width = 4, height = 4)


alpha_stats %>% 
  filter(sample_age == 1) %>% 
  wilcox_test(chao1 ~ disease_category)

alpha_stats %>% 
  filter(sample_age == 1) %>%
  ggplot(aes(x = disease_category, y = shannon, color = disease_category)) +
  geom_point(position = position_jitter(0.2), size = 3, 
             alpha = 0.5) +
  stat_summary(geom = 'bar', fun = mean, fill = NA, 
               color = 'black', width = 0.4, linewidth = 1.5) +
  stat_summary(geom = 'errorbar', 
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               color = 'black', width = 0.2, linewidth = 1.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 8)) +
  # scale_x_discrete(labels = c('Healthy', 'Acute GI')) +
  scale_color_manual(values = c( 'blue', 'red')) +
  ylab('Shannon Index') +
  theme(
    aspect.ratio = 2/1,
    legend.position = 'none',
    axis.title.x = element_blank()
  )

ggsave('plots/acute_gi/shannon_timepoint_1.png', width = 4, height = 4)

alpha_stats %>% 
  filter(sample_age == 1) %>%
  rstatix::wilcox_test(shannon ~ disease_category)



bc_dist <- read_rds('data/240711_merged/data/bc_dist_raw.Rds') %>% 
  usedist::dist_subset(., idx = alpha_stats$sampleid)
d1_samples <- alpha_stats %>% 
  filter(sample_order == 1)
d1_dist <- bc_dist %>% 
  usedist::dist_subset(., d1_samples$sampleid)

adonis2(d1_dist ~ disease_category, permutations = 9999, data = d1_samples)

d1_pcoa <- generate_pcoa(d1_dist, d1_samples)

d1_pcoa[[1]] %>% 
  ggplot(aes(x= PCo1, y = PCo2, color = disease_category)) +
  geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = c('blue', 'red'), name = 'Diagnosis') +
  xlab(glue::glue('PCo1 - {round(d1_pcoa[[2]][1],2)}%')) +
  ylab(glue::glue('PCo2 - {round(d1_pcoa[[2]][2],2)}%')) +
  theme(
    aspect.ratio = 1/1.25,
    legend.text = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold')
  )
ggsave('plots/acute_gi/d1_pcoa.png', width = 6, height = 4)

  


long_alpha_data <- alpha_stats %>%
  mutate(disease_category = case_when(disease_category == 'Healthy' ~ 'Healthy',
                                      TRUE ~ 'Acute GI')) %>% 
  mutate(horse_visit = paste0(horse, visit)) %>% 
  group_by(horse_visit) %>% 
  filter(sample_age <= 7) %>% 
  filter(n() >= 3) %>% 
  ungroup()

long_alpha_data %>% 
  group_by(disease_category) %>% 
  summarize(n = n(),
            patients = n_distinct(horse, visit))

long_alpha_data %>% 
  ggplot(aes(x = sample_age, y = chao1,  
             color = fct_relevel(disease_category, 'Healthy', after = 0))) +
  geom_line(aes(group = horse_visit), alpha = 0.2) +
  stat_summary(geom = 'line', aes(group = disease_category), fun = 'mean',
               linewidth = 1.5) +
  scale_color_manual(values = c( 'blue', 'red'),
                     labels = c('Healthy', 'Acute GI')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(0, 3000)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(1, 7), 
                     breaks = c(seq(1, 7, 1))) +
  ylab('Chao1 Index') +
  xlab('Day') +
  theme(
    aspect.ratio = 2/2.5,
    legend.text = element_text(face = 'bold')
  )

ggsave('plots/acute_gi/long_choa1.png', width = 8, height = 3)

long_alpha_data %>% 
  ggplot(aes(x = sample_age, y = shannon,              
             color = fct_relevel(disease_category, 'Healthy', after = 0))) +
  geom_line(aes(group = horse_visit), alpha = 0.3) +
  stat_summary(geom = 'line', aes(group = disease_category), fun = 'mean',
               linewidth = 1.5) +
  scale_color_manual(values = c( 'blue', 'red'),
                     labels = c('Healthy', 'Acute GI')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(0, 8)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(1, 7), 
                     breaks = c(seq(1, 7, 1))) +
  ylab('Shannon Index') +
  xlab('Day') +
  theme(
    aspect.ratio = 2/2.5,
    legend.text = element_text(face = 'bold')
  )

ggsave('plots/acute_gi/long_shannon.png', width = 8, height = 3)


# linear mixed-effects models
library(lmerTest)  
library(emmeans)

colnames(long_alpha_data)
long_alpha_data$sample_age <- as.double(long_alpha_data$sample_age)
chao1_model <- lmer(chao1 ~  sample_age + (1 | horse_visit), 
                data = long_alpha_data %>% 
                  filter(disease_category == 'Acute GI'))
anova(chao1_model)

shannon_model <- lmer(shannon ~ sample_age + (1 | horse_visit), 
                    data = long_alpha_data %>% 
                      filter(disease_category == 'Healthy'))
anova(shannon_model)


library(usedist)
library(rstatix)
source('code/data_beta_diversity.R')


bc_dist_1wk <- bc_dist %>% 
  dist_subset(., long_alpha_data$sampleid)

bc_pcoa_1wk <- generate_pcoa(bc_dist_1wk, metadata_df = long_alpha_data)

adonis2(bc_dist_1wk ~ disease_category * sample_age, permutations = 9999,
        data = long_alpha_data)


bc_pcoa_1wk[[1]] %>% 
  ggplot(aes(x = PCo1, y = PCo2, color = fct_rev(disease_category),
             alpha = sample_age)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('blue', 'red'), name = 'Diagnosis') +
  scale_alpha_continuous(breaks= c(1, 3, 5, 7), name = 'Days') +
  xlab(glue::glue('PCo1 - {round(bc_pcoa_1wk[[2]][1],2)}%')) +
  ylab(glue::glue('PCo2 - {round(bc_pcoa_1wk[[2]][2],2)}%')) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme(
    aspect.ratio = 1,
    legend.title = element_text(face = 'bold', size = 16),
    legend.text = element_text(face = 'bold',size = 14)
  )
ggsave('plots/acute_gi/acute_gi_main_dx_1wk.png', width = 5, height = 5)

t1_samples <- acute_gi_samples %>% 
  filter(sample_age == 1)

t1_samples <- long_alpha_data %>% 
  filter(sample_age == 1) %>% 
  select(sampleid, horse) %>% 
  rename(horse_id_a = 'horse')

t1_samples_b <- long_alpha_data %>% 
  rename(horse_id_b = 'horse',
         sampleid_b = 'sampleid') %>% 
  select(sampleid_b, horse_id_b)

long_metadata <- long_alpha_data %>% 
  select(sampleid, horse, visit, sample_age, sample_order, disease_category)


dist_from_d1 <- bc_dist_1wk %>% 
  as.matrix() %>% 
  as.tibble(rownames = 'sampleid') %>% 
  pivot_longer(-sampleid, names_to = 'sampleid_b', values_to = 'dist') %>% 
  filter(sampleid_b %in% t1_samples$sampleid) %>% 
  left_join(., long_metadata, by = 'sampleid') %>% 
  left_join(., t1_samples_b) %>% 
  filter(horse == horse_id_b)
  

dist_from_d1 %>% 
  group_by(horse, sample_age,  disease_category) %>% 
  summarize(dist = mean(dist)) %>% 
  ggplot(aes(x = sample_age, y = dist, 
             color = fct_relevel(disease_category, 'Healthy', after = 0))) +
  geom_line(aes(group = horse), alpha =0.2) +
  stat_summary(geom = 'line', aes(group = disease_category), fun = 'mean',
               linewidth = 1.5) +
  scale_color_manual(values = c( 'blue', 'red'),
                     labels = c('Healthy', 'Acute GI')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(0, 1)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(1,7), 
                     breaks = c(seq(1, 7, 1))) +
  ylab('Bray-Curtis Dissimilarity') +
  xlab('Day') +
  theme(
    aspect.ratio = 2/2.5,
    legend.text = element_text(face = 'bold')
  )

ggsave('plots/acute_gi/bc_dist_from_t1.png', width = 8, height = 3)

dist_from_d1_model <- lmerTest::lmer(dist ~sample_age + (1 | horse), 
                      data = dist_from_d1 %>% 
                        filter(sample_age > 1) %>% 
                        filter(disease_category == 'Acute GI'))
anova(dist_from_d1_model)


library(mikropml)

adult_metadata <- acute_gi_samples %>% 
  filter(sample_age == 1) %>% 
  mutate(age_of_collection = as.double(age_of_collection)/365) %>%
  filter(age_of_collection >= 4) %>% 
  select(sampleid, disease_category) %>% 
  mutate(disease_category = case_when(disease_category == 'Acute GI' ~ 'Acute_GI',
                                      TRUE ~ disease_category)) 

g_names <- taxonomy %>% 
  pivot_wider(names_from = 'level', values_from = 'taxon') %>% 
  mutate(Genus = paste(Phylum, Class, Order, Family, Genus, sep = '_')) %>% 
  select(featureid, Genus) %>% 
  rename(taxon = 'Genus') 

g_table <- table %>% 
  select(featureid, adult_metadata$sampleid) %>% 
  pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
  left_join(., g_names) %>% 
  group_by(taxon, sampleid) %>% 
  summarize(count = sum(count), .groups = 'drop') %>% 
  group_by(sampleid) %>% 
  mutate(rel_abund = count / sum(count)) %>% 
  ungroup() %>% 
  select(-count) %>% 
  pivot_wider(names_from = 'taxon', values_from = 'rel_abund') 


model_data <- g_table %>% 
  filter(sampleid %in%  adult_metadata$sampleid) %>% 
  left_join(., adult_metadata) %>% 
  column_to_rownames(var = 'sampleid')

model_data_processed <- preprocess_data(dataset = model_data, outcome_colname = 'disease_category')

as.data.frame(model_data_processed$removed_feats) %>% 
  write_csv(glue::glue('data/acute_gi/rf_model/pre_processing/day1_removed_features.csv'))
as.data.frame(model_data_processed$dat_transformed) %>% 
  write_csv(glue::glue('data/acute_gi/rf_model/pre_processing/day1_transformed_data.csv'))

d1_rf_model_genus_rel_abund <- run_ml(dataset = model_data_processed$dat_transformed,
                                      method = 'rf',
                                      outcome_colname = "disease_category",
                                      seed = 1851, 
                                      find_feature_importance = T
)

d1_rf_model_genus_rel_abund$trained_model %>% 
  saveRDS(glue::glue('data/acute_gi/rf_model/models/d1_genus_rel_abund_rf_model.Rds'))

d1_rf_model_genus_rel_abund$performance

# model_list <- list.files('data/acute_gi/rf_model/models/')
d1_rf_model_genus_rel_abund %>% 
  pluck("performance") %>% 
  write.csv('data/acute_gi/rf_model/performance/d1_performance_summary.csv', row.names = F)

features <- d1_rf_model_genus_rel_abund$feature_importance %>% 
  as_tibble() %>% 
  arrange(-perf_metric_diff) 

features %>% 
  write.csv(glue::glue('data/acute_gi/rf_model/feature_importance/rf_model_genus_rel_abund_d1_features.csv'), row.names = F)

plot <- features %>% 
  filter(pvalue < 0.05) %>% 
  ggplot(aes(x = fct_reorder(feat, perf_metric_diff), y = perf_metric_diff)) +
  geom_col() +
  coord_flip() +
  ggprism::theme_prism() +
  scale_y_continuous(expand = expansion(mult = c(0,0.25))) +
  theme(
    axis.title.y = element_blank(),
  )

ggsave('data/acute_gi/rf_model/feature_importance/rf_model_genus_rel_abund_d1_features.png', plot = plot,
       width = 12)

d1_rf_model_genus_rel_abund$performance

sensspec <- calc_model_sensspec(
  trained_model = d1_rf_model_genus_rel_abund$trained_model,
  test_data = d1_rf_model_genus_rel_abund$test_data, 
  outcome_colname = 'disease_category')

mean_roc <- calc_mean_roc(sensspec)

mean_roc %>% 
  ggplot(aes(x = 1 -specificity, y = mean_sensitivity)) +
  geom_line(linewidth = 1, color = 'red')  +
  geom_segment(x = 0, y = 0, xend = 1, yend = 1, linetype = 2,
               inherit.aes = F, linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), 
                     limits = c(0,1),
                     name = 'Sensitivity') +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.01)), limits = c(0,1),
                     name = '1 - Specificity')  +
  theme(aspect.ratio = 1)
ggsave(glue::glue('data/acute_gi/rf_model/performance/rf_model_genus_rel_abund_d1_roc.png'), width = 5, height = 4)


d1_features <- features

d1_features_enriched <- g_table %>% 
  filter(sampleid %in% adult_metadata$sampleid) %>% 
  pivot_longer(-sampleid, names_to = 'feat') %>% 
  left_join(., adult_metadata) %>% 
  filter(feat %in% d1_features$feat) %>% 
  group_by(feat, disease_category) %>% 
  summarise(mean = mean(value),.groups = 'drop') %>% 
  pivot_wider(names_from = 'disease_category', 
              values_from = 'mean') %>% 
  mutate(enriched = ifelse(Acute_GI > Healthy, 'Acute_GI', 'Healthy')) %>% 
  select(feat, enriched)

shortened_feature_name <- c('Escherichia-Shigella', 'Bacteroides', 'Oscillospirales UCG-010',
                            'Lachnospiraceae UCG-006', 'Lachnospiraceae UCG-009',
                            'Butyricicoccaceae UCG-009', 'Pseudobutyrivibrio',
                            'Bacteroidales', 'Uncultured Bacteroidales', 'Oscillospiraceae',
                            'Streptococcus', 'Clostridium methylpentosum Gp.', 
                            'Eubacterium nodatum Gp.', 'Roseburia', 'Endomicrobium'
                            )
                              


d1_features %>% 
  left_join(., d1_features_enriched) %>% 
  mutate(enriched = factor(enriched, levels = c('Healthy', 'Acute_GI'))) %>% 
  top_n(n = 15, wt = perf_metric_diff) %>%
  arrange(-perf_metric_diff) %>% 
  cbind(., shortened_feature_name) %>%
  # select(feat, shortened_feature_name)
  ggplot(aes(x = fct_reorder(shortened_feature_name, perf_metric_diff), y = perf_metric_diff,
             fill = (enriched))) +
  geom_col(color = NA) +
  ggprism::theme_prism()  +
  scale_y_continuous(expand = expansion(mult = c(0,0.25))) +
  scale_fill_manual(values = c('blue', 'red'), labels = c('Healthy', 'Acute GI'),
                    name = 'Enriched') +
  ylab('Importance') +
  theme(
    aspect.ratio = 1,
    axis.title.y = element_blank(),
    legend.text = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold', size = 14),
    legend.position = c(0.75, 0.2)
  ) +
  coord_flip() 

ggsave('data/acute_gi/rf_model/feature_importance/d1_top_15_features.png', height = 3.5, width = 7)

d1_features %>% 
  filter(pvalue < 0.05) %>% 
  arrange(-perf_metric_diff) %>% 
  print(n = 100)

d1_features %>% 
  left_join(., d1_features_enriched) %>% 
  mutate(enriched = factor(enriched, levels = c('Healthy', 'Acute_GI'))) %>% 
  ggplot(aes(x = perf_metric_diff, y = -log10(pvalue), color = enriched)) +
  geom_point(size = 3) +
  geom_vline(aes(xintercept = 0), linetype = 2, linewidth = 1) +
  geom_hline(aes(yintercept = -log10(0.05)), linetype = 2, linewidth = 1) +
  scale_color_manual(values = c('Healthy' = 'blue', 'Acute_GI' =  'red','NA' =  'gray'),
                     breaks = c('Healthy', 'Acute_GI'),
                     labels = c('Healthy', 'Acute GI'),
                     name = 'Diagnosis') +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  xlab('Importance') +
  ylab('-Log<sub>10</sub>(<i>p</i> value)') +
  theme(
    aspect.ratio = 1/1.25,
    axis.title.y = ggtext::element_markdown(),
    legend.text = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold'),
    legend.position = c(0.8, 0.15)
  )
ggsave('data/acute_gi/rf_model/feature_importance/d1_features.png', width = 5, height = 4)


### Subgroups
acute_gi_groups <- alpha_stats %>% 
  mutate(acute_gi_subgroup = case_match(acute_gi_subgroup,
                                        NA ~ 'Healthy',
                                        .default = acute_gi_subgroup),
         colic_inflammatory = case_match(colic_inflammatory,
                                         NA ~ 'Healthy',
                                         .default = colic_inflammatory),
         inflam_group = paste(colic_inflammatory, acute_gi_subgroup, sep = '_'))

acute_gi_groups %>% 
  # filter(sample_age == 1) %>% 
  group_by(colic_inflammatory, acute_gi_subgroup) %>% 
  summarise(samples = n(), 
            subject = n_distinct(horse)) 

acute_gi_groups %>% 
  group_by(horse, visit) %>% 
  filter(n() >= 3) %>% 
  group_by(colic_inflammatory, acute_gi_subgroup) %>% 
  summarise(samples = n(), 
            subject = n_distinct(horse)) 


unique(acute_gi_groups$inflam_group)

pal <- c(
  "no_Other" = '#d2b4de',
  "no_Non-specific" = '#4a235a', 
  "no_Mechanical" = '#884ea0',
  "yes_Mechanical" = '#145a32',
  "yes_Infectious" = '#82e0aa',
  "yes_Other" = '#27ae60',     
  "Healthy_Healthy" = '#2980b9'
)

acute_gi_groups %>% 
  group_by(colic_inflammatory, inflam_group) %>% 
  summarize(n = n_distinct(horse, visit)) %>% 
  ggplot(aes(x = fct_reorder(inflam_group, n), 
             y = n,
             fill = inflam_group)) +
  geom_col() +
  facet_grid(colic_inflammatory~., scales = 'free', space = 'free') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = pal) +
  scale_x_discrete(breaks = c("Healthy_Healthy", "no_Non-specific",
                              "no_Mechanical", "no_Other",
                              "yes_Mechanical", "yes_Other",
                              "yes_Infectious"
                              ),
    labels = c('Healthy', 'Idiopathic', 'Mechanical',
                              'Other', 'Mechanical', 'Other', 'Infectious')
    ) +
  coord_flip() +
  ylab('# of Subjects') +
  theme(
    strip.text = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    axis.title.y = element_blank()
  )
ggsave('plots/acute_gi/subgroup_patients.png', width = 2.75, height = 4.5)

acute_gi_groups %>% 
  filter(sample_age == 1) %>% 
  ggplot(aes(x = fct_rev(fct_infreq(inflam_group, chao1,)),
             y = chao1, color = inflam_group)) +
  geom_point(position = position_jitter(), size = 2) +
  stat_summary(geom = 'bar', fun = 'mean', fill = NA, color = 'black', linewidth = 1) +
  stat_summary(geom = 'errorbar', linewidth = 1, color = 'black', width = 0.5,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x)) +
  facet_grid(colic_inflammatory~., scales = 'free', space = 'free') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)),
                     breaks = c(0, 1500, 3000)) +
  scale_color_manual(values = pal) +
  scale_x_discrete(breaks = c("Healthy_Healthy", "no_Non-specific",
                              "no_Mechanical", "no_Other",
                              "yes_Mechanical", "yes_Other",
                              "yes_Infectious"
  ),
  labels = c('Healthy', 'Idiopathic', 'Mechanical',
             'Other', 'Mechanical', 'Other', 'Infectious')) +
  coord_flip() +
  ylab('Chao1 Index') +
  theme(
    strip.text = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    axis.title.y = element_blank()
  )
ggsave('plots/acute_gi/subgroup_chao1.png', width = 3.5, height = 4.5)

acute_gi_groups %>% 
  group_by(acute_gi_subgroup) %>% 
  filter(sample_age == 1) %>% 
  filter(n() >= 3) %>% 
  ungroup() %>% 
  kruskal_test(chao1 ~ inflam_group)
acute_gi_groups %>% 
  group_by(acute_gi_subgroup) %>% 
  filter(sample_age == 1) %>% 
  filter(n() >= 3) %>% 
  ungroup() %>% 
  pairwise_wilcox_test(chao1 ~ inflam_group, p.adjust.method = 'BH') %>% 
  write_csv('plots/acute_gi/shannon_pairwise_comparison_wilcoxon.csv')

source('code/data_beta_diversity.R')
bc_dist <- read_rds('data/240711_merged/data/bc_dist_raw.Rds')
acute_gi_groups_dist <- bc_dist %>% 
  usedist::dist_subset(., acute_gi_groups %>% filter(sample_age == 1) %>% pull(sampleid))
acute_gi_groups_pcoa <- generate_pcoa(dist = acute_gi_groups_dist, metadata_df = acute_gi_groups)

acute_gi_groups_pcoa[[1]] %>% 
  ggplot(aes(x = PCo1, y = PCo2, 
             color = factor(inflam_group, levels = c("Healthy_Healthy", "no_Mechanical",
                                                     "no_Non-specific", "no_Other",
                                                     "yes_Mechanical", "yes_Infectious",
                                                     "yes_Other")))) +
  geom_point(size = 3)  +

  scale_color_manual(values = pal, name = 'Colic\nSubgroup', 
                     breaks = c("Healthy_Healthy", "no_Non-specific",
                                "no_Mechanical", "no_Other",
                                "yes_Mechanical", "yes_Other",
                                "yes_Infectious"
                     ),
                     labels = c('Healthy', 'Idiopathic', 'Mechanical',
                                'Other', 'Mechanical', 'Other', 'Infectious')) +
  xlab(glue::glue('PCo1 - {round(acute_gi_groups_pcoa[[2]][1], 2)}%')) +
  ylab(glue::glue('PCo2 - {round(acute_gi_groups_pcoa[[2]][2], 2)}%')) +
  theme(
    aspect.ratio = 2/1.5,
    legend.text = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold', size = 14, hjust = 0.5)
  )

ggsave('plots/acute_gi/subgroup_pcoa_d1.png', width = 6, height = 4.5)
 

acute_gi_subgroup_tmp <- acute_gi_groups %>% 
  group_by(inflam_group) %>% 
  filter(sample_age == 1) %>% 
  filter(n() >=3)

acute_gi_groups_dist <- bc_dist %>% 
  usedist::dist_subset(., acute_gi_subgroup_tmp$sampleid)

adonis2(acute_gi_groups_dist ~ inflam_group, permutations = 9999,
        data = acute_gi_subgroup_tmp)

pairwiseAdonis::pairwise.adonis2(acute_gi_groups_dist ~ inflam_group, nperm = 9999, 
                                 data = acute_gi_subgroup_tmp)


subgroup_long <- acute_gi_samples %>% 
  left_join(., alpha_stats) %>% 
  group_by(horse, visit) %>% 
  filter(sample_age <= 5) %>%
  filter(n() >= 3) %>% 
  mutate(long_group = case_when(colic_inflammatory == 'yes' ~ 'Inflammatory',
                                colic_inflammatory == 'no' & 
                                  acute_gi_subgroup == 'Non-specific' ~ 'Idiopathic',
                                colic_inflammatory == 'no' & 
                                  acute_gi_subgroup == 'Mechanical' ~ 'Mechanical',
                                disease_category == 'Healthy' ~ 'Healthy')) %>% 
  ungroup() %>% 
  filter(!is.na(long_group)) %>% 
  mutate(long_group = factor(long_group, levels = c('Healthy', 'Idiopathic', 'Mechanical', 'Inflammatory')))

subgroup_long %>% 
  ggplot(aes(x = sample_age, y = shannon, group = horse, color = long_group)) +
  geom_line(alpha = 0.3) +
  stat_summary(geom = 'line', fun = 'mean', linewidth = 1.5,
               aes(group = long_group)) +
  scale_color_identity(name = 'Colic\nSubgroup', guide = 'legend') +
  scale_color_manual(values = c('#2980b9', '#4a235a', '#884ea0', '#27ae60'),
                     name = 'Colic\nSubgroup')+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 8)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(1,5), breaks = seq(1,7,1)) +
  xlab('Day') +
  ylab('Shannon Index') +
  theme(
    aspect.ratio = 1.25/2,
    legend.text = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold', hjust = 0.5),
    # legend.position = 'none'
  )

ggsave('plots/acute_gi/subgroup_long_chao1.png', width = 7, height = 4.5)


library(lmerTest)

subgroup_long_lm<- lmer(chao1 ~ long_group * sample_age + (1|horse), 
                      data = subgroup_long)
anova(subgroup_long_lm)

subgroup_long_lm <- lmer(shannon ~ long_group * sample_age + (1|horse), 
                         data = subgroup_long)
anova(subgroup_long_lm)

sugroup_long_main <- subgroup_long %>% ungroup() 
long_dist <- bc_dist %>% 
  usedist::dist_subset(., idx = subgroup_long$sampleid)

long_pcoa <- generate_pcoa(long_dist, metadata_df = subgroup_long)


pal = c('Inflammatory' = '#27ae60',
        'Idiopathic' = '#4a235a',
        'Mechanical' = '#884ea0',
        'Healthy' = '#2980b9')
long_pcoa[[1]] %>% 
  ggplot(aes(x = PCo1, y = PCo2, color = long_group, alpha = sample_age,
             )) +
  geom_point(size = 6, color = 'gray', alpha = 0.5,
             data = long_pcoa[[1]] %>% filter(!long_group %in% c('Healthy', 'Idiopathic'))) +
  geom_point(size = 6,
             data = long_pcoa[[1]] %>% filter(long_group %in% c('Healthy'))) +
  geom_point(size = 6,
             data = long_pcoa[[1]] %>% filter(long_group %in% c('Idiopathic'))) +
  scale_color_manual(values = pal, name = 'Colic\nSubgroup') +
  scale_alpha_continuous(name = 'Sample\nAge', breaks = c(1, 3, 5, 7), range = c(0.25, 1)) +
  xlab(glue::glue('PCo1 - {round(long_pcoa[[2]][1], 2)}%')) +
  ylab(glue::glue('PCo2 - {round(long_pcoa[[2]][2], 2)}%')) +
  theme(
    aspect.ratio = 1,
    legend.title = element_text(face = 'bold', hjust = 0.5),
    legend.text = element_text(face = 'bold', size = 14)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 4)) ,
    alpha = guide_legend(order = 2, override.aes = list(size = 4)) ,
    
  )

ggsave('plots/acute_gi/healthy_Mechanical_pcoa.png', width = 6, height = 4)

adonis2(long_dist ~ long_group * sample_age, permutations = 9999,
        data = long_pcoa[[1]])
pairwiseAdonis::pairwise.adonis2(long_dist ~ long_group * sample_age, nperm = 9999,
                                 data = sugroup_long_main)


# BC Distance from Day 1 sample
t1_samples <- subgroup_long %>% 
  filter(sample_age == 1) %>% 
  select(sampleid, horse) %>% 
  rename(horse_id_a = 'horse')

t1_samples_b <- subgroup_long %>% 
  rename(horse_id_b = 'horse',
         sampleid_b = 'sampleid') %>% 
  select(sampleid_b, horse_id_b)

long_metadata <- subgroup_long %>% 
  select(sampleid, horse, visit, sample_age, sample_order, long_group)

dist_from_d1 <- long_dist %>% 
  as.matrix() %>% 
  as.tibble(rownames = 'sampleid') %>% 
  pivot_longer(-sampleid, names_to = 'sampleid_b', values_to = 'dist') %>% 
  filter(sampleid_b %in% t1_samples$sampleid) %>% 
  left_join(., long_metadata, by = 'sampleid') %>% 
  left_join(., t1_samples_b) %>% 
  filter(horse == horse_id_b)

dist_from_d1 %>% 
  group_by(long_group, sample_age) %>% 
  filter(long_group %in% c('Healthy', 'Inflammatory')) %>%
  ggplot(aes(x = sample_age, y = dist, group = horse, color = long_group)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = long_group), linewidth = 1.5, geom = 'line') +
  scale_color_manual(values = pal, name = 'Subgroup') +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(1,7,1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  xlab('Day') +
  ylab('Bray-Curtis Dissimilarity') +
  theme(
    aspect.ratio = 1,
    legend.title = element_text(face = 'bold', hjust = 0.5),
    legend.text = element_text(face = 'bold', size = 14),
    legend.position = 'none'
  ) +
  guides(
    color = guide_legend(order = 1)  
  )
ggsave('plots/acute_gi/dist_from_d1.png', width = 4, height = 4)


