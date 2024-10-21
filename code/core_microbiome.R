library(tidyverse)
library(microbiome)
library(RColorBrewer)
library(rstatix)
library(ggprism)
library(phyloseq)
library(microViz)
library(randomForest)
library(lmerTest)

theme_set(theme_prism())

source('code/load_data.R')


healthy_adult_metadata <- metadata %>% 
  mutate(age_of_collection = as.double(age_of_collection) / 365,
         month = case_when(country == 'ASTL' & month_norm <= 6 ~ month_norm + 6,
                           country == 'ASTL' & month_norm > 6 ~ month_norm - 6,
                           TRUE ~ month_norm)) %>% 
  filter(disease_category == 'Healthy' &
           age_of_collection >= 4) %>% 
  filter(sampleid != 'JW1')  

alpha_stats <- read_rds('data/240711_merged/data/alpha.stats.Rds') %>% 
  filter(sampleid %in% healthy_adult_metadata$sampleid) %>% 
  left_join(healthy_adult_metadata)


anova(lmer(chao1 ~ as.double(sample_age) + (1 | horse), data = alpha_stats))
anova(lmer(shannon ~ as.double(sample_age) + (1 | horse), data = alpha_stats))


chao1_data <- alpha_stats %>% 
  mutate(age = as.double(age_of_collection)) %>% 
  select(sampleid, chao1, sex, breed_type, age, 
         country,month_norm) %>% 
  drop_na() %>% 
  column_to_rownames(var = 'sampleid')

shannon_data %>% 
  rstatix::cor_test(shannon, age, method = 'spearman' )

chao1_rf_model <- randomForest(chao1 ~ ., data = chao1_data, proximity = T, 
                               importance = T, ntree = 1000)
plot(chao1_rf_model)
importance(chao1_rf_model)

chao1_importance <- importance(chao1_rf_model) %>% 
  as_tibble(rownames = 'category') %>% 
  select(category, `%IncMSE`) %>% 
  rename(IncMSE_chao1 = `%IncMSE`) 

chao1_importance %>% 
  mutate(category = str_replace_all(category, '_', ' '),
         category = str_to_title(category)) %>% 
  ggplot(aes(x = fct_reorder(category, IncMSE_chao1), y = IncMSE_chao1)) +
  geom_col(fill = 'blue')+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, 31),
                     breaks = c(0, 10, 20, 30)) +
  ggprism::theme_prism() +
  scale_x_discrete(labels = c( 'Sex', 'Month','Breed Type','Age', 'Country')) +
  ylab('%IncMSE') +
  theme(
    aspect.ratio = 3/1.25,
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  coord_flip() 

ggsave('plots/fig2_ideas/rf_importance_chao1.png', width = 3, height = 3)

shannon_data <- alpha_stats %>% 
  mutate(age = as.double(age_of_collection)) %>% 
  select(sampleid, shannon, sex, breed_type, age, month_norm, 
         country) %>% 
  drop_na() %>%   
  column_to_rownames(var = 'sampleid')

shannon_rf_model <- randomForest(shannon ~ ., data = shannon_data, proximity = T, 
                                 importance = T, ntree = 1000)
plot(shannon_rf_model)
importance(shannon_rf_model)
shannon_importance <- importance(shannon_rf_model) %>% 
  as_tibble(rownames = 'category') %>% 
  select(category, `%IncMSE`) %>% 
  rename(IncMSE_shannon = `%IncMSE`) 

shannon_importance %>% 
  mutate(category = factor(category, 
                           levels = rev(c('country', 'age', 'breed_type', 'month_norm', 'sex')))) %>% 
  ggplot(aes(x = category,
             y = IncMSE_shannon)) +
  geom_col(fill = 'blue')  +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, 32),
                     breaks = c(0, 10, 20, 30)) +
  ggprism::theme_prism() +
  ylab('%IncMSE') +
  theme(
    aspect.ratio = 3/1.25,
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  coord_flip() 

ggsave('plots/fig2_ideas/rf_importance_shannon.png', width = 3, height = 3)


chao1_importance %>% 
  left_join(., shannon_importance) %>%
  # arrange(-IncMSE_chao1)
  write_csv('plots/fig2_ideas/alpha_stat_rd_incMSE.csv')

choa1_data %>% 
  pairwise_wilcox_test(chao1 ~ country)
choa1_data %>% 
  ggplot(aes(x = country, y = chao1, color = country)) +
  geom_point(size = 3, position = position_jitter()) +
  stat_summary(geom = 'bar', fun = 'mean', fill = NA, 
               color = 'black', linewidth = 1) +
  stat_summary(geom = 'errorbar',
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x), 
               color = 'black', linewidth = 1, width = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)),
                     breaks = c(0, 1000, 2000)
  ) +
  ylab('Chao1 Index') +
  scale_color_manual(values = rev(c('blue', 'red', 'darkgreen'))) +
  theme(
    aspect.ratio = 3/1.25,
    legend.position = 'none',
    axis.title.y = element_blank()
  ) +
  coord_flip()

ggsave('plots/fig2_ideas/chao1_country.png', width = 4, height = 3.5)

choa1_data %>% 
  kruskal_test(chao1 ~ country)
choa1_data %>% 
  pairwise_wilcox_test(chao1 ~ country, p.adjust.method = 'BH')



source('code/data_beta_diversity.R')
bc_dist <- read_rds('data/240711_merged/data/bc_dist_raw.Rds')


healthy_bc_dist <- bc_dist %>% 
  usedist::dist_subset(., healthy_adult_metadata$sampleid)
healthy_pcoa <- generate_pcoa(healthy_bc_dist, healthy_adult_metadata)

# adonis2(healthy_bc_dist ~ sex +  age_of_collection  + month_norm + country + breed_type, permutations = 9999,
#         data = healthy_adult_metadata)
# tmp <- .Last.value %>%
#   as_tibble(rownames = 'effect') 

healthy_permanova_res <- read_csv('plots/fig2_ideas/healthy_permanova_res.csv')

healthy_permanova_res %>% 
  drop_na() %>% 
  # rename(effect = `effect...3`)%>%
  mutate(effect = case_when(effect == 'age_of_collection' ~ 'age', TRUE ~ effect)) %>% 
  ggplot(aes(x = factor(effect, levels = rev(c('country', 'breed_type', 'age',
                                               'sex', 'month_norm'))), y = R2)) +
  geom_col(fill = 'blue') +
  coord_flip() +
  ylab('R^2^') +
  scale_y_continuous(limits = c(0, 0.1),
                     expand = expansion(mult = c(0, 0.05)),
                     breaks = c(0,0.05, 0.1),
                     labels = scales::percent) +
  # scale_x_discrete(labels = c('Month', 'Age', 'Sex', 'Country', 'Breed Type')) +
  theme(
    aspect.ratio = 3/1.25,
    axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    axis.title.x = ggtext::element_markdown(),
    axis.ticks.y = element_blank()
  )
ggsave('plots/fig2_ideas/Healthy_permanova_res.png', width = 3, height = 3)



healthy_pcoa[[1]] %>% 
  ggplot(aes(x = PCo1, y = PCo2, color = country)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('darkgreen', 'red', 'blue'), 
                     name = 'Country') +
  xlab(glue::glue('PCo1 - {round(healthy_pcoa[[2]][1],2)}%')) +
  ylab(glue::glue('PCo2 - {round(healthy_pcoa[[2]][2],2)}%')) +
  theme(
    aspect.ratio = 1,
    legend.text = element_text(face = 'bold', size = 12),
    legend.title = element_text(face = 'bold', size = 12)
  ) 
ggsave('plots/fig2_ideas/country_pcoa.png', width = 4.5, height = 3.5)


breed_colors <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", 
                           "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
                           
healthy_pcoa[[1]] %>% 
  ggplot(aes(x = PCo1, y = PCo2, color = breed_type)) +
  geom_point(size = 3) +
  scale_color_manual(values = rev(breed_colors),
                     breaks = c('Draught', 'Cob', 'Pony',
                                'Warmblood', 'Hot-blooded', 'Other',
                                'Quarter Horse', 'Donkey'),
                     name = 'Breed Type') +
  xlab(glue::glue('PCo1 - {round(healthy_pcoa[[2]][1],2)}%')) +
  ylab(glue::glue('PCo2 - {round(healthy_pcoa[[2]][2],2)}%')) +
  theme(
    aspect.ratio = 1,
    legend.text = element_text(face = 'bold', size = 12),
    legend.title = element_text(face = 'bold', size = 12)
  ) 

ggsave('plots/fig2_ideas/breed_type_pcoa.png', width = 4.5, height = 3.5)


healthy_adult_metadata %>% 
  select(horse, breed_type) %>% 
  distin
metadata %>% 
  filter(disease_category == 'Healthy') %>% 
  filter(age_of_collection > 365*4) %>% 
  # select(horse, breed_type) %>% 
  distinct()


## Main Figure 3

physeq <- read_rds('data/240711_merged/physeq.Rds')

# Following tutorial from https://microbiome.github.io/tutorials/Core.html
healthy_all <- ps_filter(physeq, disease_category == 'Healthy', age_of_collection >= 365 * 4)
healthy_family <- tax_glom(healthy_all, taxrank = 'Family')
healthy_rel <- microbiome::transform(healthy_family, "compositional")

prevalences <- seq(.05, 1, .05)
detections <- seq(0.01, 0.1, 0.01)

core <- plot_core(healthy_rel,
                  plot.type = "heatmap", 
                  prevalences = prevalences, 
                  detections = detections, 
                  min.prevalence = prevalence(healthy_rel, sort = TRUE)[100]) 

keep_taxonomy <- core$data %>%
  distinct(Taxa)

taxa_label<- taxonomy %>% 
  filter(level %in% c("Phylum",'Order', "Family")) %>% 
  rename(Taxa = "featureid") %>% 
  filter(Taxa %in% keep_taxonomy$Taxa) %>% 
  pivot_wider(values_from = 'taxon', names_from = 'level') %>% 
  mutate(Family = str_replace_all(Family, '_', ' '),
         Family = str_replace(Family, '\\[', ''),
         Family = str_replace(Family, '\\]', ''),
         Family = str_replace(Family, 'group', 'Gp.'),
         Family = case_match(Family,
                             'UCG-010' ~ 'Oscillospirales UCG-010',
                             .default = Family)
         ) %>% 
  select(-Order)

core$data <- core$data %>% 
  left_join(., taxa_label, by = "Taxa") %>% 
  mutate(Taxa = Family) %>% 
  select(Taxa, DetectionThreshold, Prevalence) %>% 
  arrange(DetectionThreshold, Prevalence)

core$data %>% 
  mutate(DetectionThreshold = as.double(DetectionThreshold) * 1) %>%
  ggplot(aes(x = as.factor(DetectionThreshold), y = fct_reorder(Taxa, Prevalence, .fun = 'mean'), 
             fill = Prevalence, color = Prevalence)) +
  geom_tile() +
  scale_fill_gradientn(colors= rev(brewer.pal(5, "Spectral")),
                       guide = guide_colorbar(ticks = F),
                       labels = scales::percent,
                       limits = c(0,1)) +
  scale_color_gradientn(colors= rev(brewer.pal(5, "Spectral")),
                       guide = guide_colorbar(ticks = F),
                       labels = scales::percent,
                       limits = c(0,1)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.text = element_text(face = "bold", color = "black", size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face = "bold", color = "black", size = 16),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.title = element_text(face = "bold", hjust = 0.5, color = "black", size = 16),
        legend.position = 'top',
        legend.title.position = 'top',
        legend.ticks = element_blank()) +
  labs(x = "Relative Abundance (%)")
ggsave('plots/fig2_ideas/core_microbiome.png', width = 10, height = 7)

order <- c(
  'WCHB1-41', 'Lachnospiraceae','Rikenellaceae',
  'Oscillospiraceae', 'Prevotellaceae', 'p-251-o5',
  'F082', 'Christensenellaceae', 'Anaerovoracaceae',
  'Oscillospirales UCG-010', 'Ruminococcaceae', 
  'Eubacterium coprostanoligenes Gp.',
  'Spirochaetaceae','Erysipelotrichaceae', 'Fibrobacteraceae', 
  'Bacteroidales UCG-001', 'RF39', 'Bacteroidales RF16 Gp.',
  'Hungateiclostridiaceae', 'Clostridia UCG-014', 'Acidaminococcaceae'
)


f_names <- taxonomy %>% 
  filter(level %in% c('Phylum', 'Order', 'Family')) %>% 
  pivot_wider(names_from = 'level', values_from = 'taxon') %>% 
  unite('Family', Phylum:Family, sep = ' ')

f_table <- table %>% 
  select(featureid, healthy_adult_metadata$sampleid) %>% 
  pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
  left_join(., f_names) %>% 
  group_by(sampleid, Family) %>% 
  summarize(count = sum(count)) %>% 
  group_by(sampleid) %>% 
  mutate(rel_abund = count / sum(count)) %>% 
  separate(Family, into = c('Phylum', 'Order', 'Family'), sep = " ")


core_family_rel_abund_table <- f_table %>%
  mutate(Family = str_replace_all(Family, '_', ' '),
         Family = str_replace(Family, '\\[', ''),
         Family = str_replace(Family, '\\]', ''),
         Family = str_replace(Family, 'group', 'Gp.'),
         Family = case_match(Family,
                             'UCG-010' ~ 'Oscillospirales UCG-010',
                             .default = Family)
  ) %>% 
  filter(sampleid %in% healthy_adult_metadata$sampleid &
           Family %in% order) %>% 
  left_join(., healthy_adult_metadata) %>% 
  select(sampleid, Phylum, Family, rel_abund, age_of_collection, breed_type, horse, visit, 
           sex, country, month_norm)


family_effect_size <- function(family){
  tmp_family_lm_res <- core_family_rel_abund_table %>% 
    filter(Family == family) %>% 
    lmer(rel_abund ~ breed_type + sex + country + month_norm + 
           age_of_collection + (1|horse),  data = .)
  
  effect_size <- effectsize::eta_squared(tmp_family_lm_res,
                                         partial = T) %>% 
    as_tibble()  %>% 
    select(Parameter, Eta2_partial)

  anova_res <- anova(tmp_family_lm_res) %>% 
    as_tibble(rownames = 'Parameter') %>% 
    select(Parameter, `Pr(>F)`) %>% 
    rename('p' = `Pr(>F)`) %>% 
    add_significance(p.col = 'p',
                     cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                     symbols = c('***', '**', '*', ''))
  
  effect_size %>% 
    left_join(., anova_res) %>% 
    mutate(family = family)
  
}

summary <- map_dfr(order, family_effect_size) 

summary %>% arrange(-Eta2_partial)

summary %>% 
  ggplot(aes(x = factor(family, levels = rev(order)), 
             y = fct_reorder(Parameter, -Eta2_partial))) +
  geom_tile(aes(fill = Eta2_partial, color = Eta2_partial), width = 1, height = 1) +
  geom_text(aes(label = p.signif), color = 'white', size = 8,
            position = position_nudge(x = -0.25)) +
  scale_fill_gradientn(colors= rev(brewer.pal(5, "Spectral")),
                       guide = guide_colorbar(ticks = F),
                       limits = c(0, 0.5), breaks = seq(0, 0.50, 0.25),
                       name = 'Partial Eta Squared') +
  scale_color_gradientn(colors= rev(brewer.pal(5, "Spectral")),
                       guide = guide_colorbar(ticks = F),
                       limits = c(0, 0.5), breaks = seq(0, 0.50, 0.25),
                       name = 'Partial Eta Squared') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0),
                   labels = c('Country', 'Breed Type', 'Age',
                              'Month', 'Sex')) +
  theme(
    aspect.ratio = 2/1,
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = 'top',
    legend.title = element_text(face = 'bold', hjust = 0.5),
    legend.text = element_text(face = 'bold'),
    legend.title.position = 'top',
    legend.ticks = element_blank()
  ) +
  coord_flip() 

ggsave('plots/fig2_ideas/core_microbiome_effect_size.png', height = 7, width = 6)


sig_country_families <- summary %>% 
  filter(Parameter == 'country') %>% 
  filter(p < 0.05) 

core_family_rel_abund_table %>% 
  filter(Family %in% sig_country_families$family) %>% 
  ggplot(aes(x = fct_rev(fct_reorder(Family, rel_abund, .fun = 'mean')), 
             y = rel_abund, fill = country)) +
  stat_summary(geom = 'bar', fun = 'mean',size = 1,
               position = position_dodge(0.7), color = 'black',
               width = 0.7) +
  stat_summary(geom = 'errorbar', aes(group = country), width = 0.4, size = 1,
               fun.min = function(x) mean(x),
               fun.max = function(x) mean(x) + sd(x), position = position_dodge(0.7)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), 
                     labels = scales::percent) +
  
  scale_fill_manual(values = c('darkgreen', 'red', 'blue'),
                    breaks = c('ASTL', 'UK', 'USA'),
                    labels = c('Australia', 'United Kingdom', 'United States')) +
  ylab('Relative Abundance') +
  theme(
    aspect.ratio = 1/2.5,
    axis.text.x = element_text(angle = 45, hjust = 1, vjust =1),
    axis.title.x = element_blank(),
    legend.position = c(0.8, 0.8),
    legend.text = element_text(face = 'bold', size=14),
  )

ggsave('plots/fig2_ideas/country_family.png', width = 10, height = 5)

core_family_rel_abund_table %>% 
  filter(Family %in% sig_country_families$family) %>% 
  group_by(Family) %>% 
  pairwise_wilcox_test(rel_abund ~ country, p.adjust.method = 'BH') %>% 
  write_csv('plots/fig2_ideas/country_sig_differences.csv')



## Supplementary Figure 5
core_family_rel_abund_table %>% 
  group_by(Family) %>% 
  mutate(prev = cume_dist(-rel_abund)) %>% 
  ggplot(aes(x = fct_reorder(Family, rel_abund, .fun = 'mean'), 
             y = rel_abund, color = prev)) +
  geom_point(position = position_jitter(0.2)) +
  stat_summary(geom = 'bar', fun = 'mean', fill = NA, color = 'black', 
               size = 1, width = 0.7) +
  stat_summary(geom = 'errorbar', size = 1, width = 0.5,
               fun.min = function(x) mean(x),
               fun.max = function(x) mean(x) + sd(x)) +
  coord_flip()+
  scale_color_gradientn(colors= rev(brewer.pal(5, "Spectral")),
                        guide = guide_colorbar(ticks = F),
                        labels = scales::percent,
                        limits = c(0,1), name = 'Prevalence') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = scales::percent) +
  labs(y = 'Relative Abundance',
       x = NULL) +
  theme(
    aspect.ratio = 2.5/1,
    legend.ticks = element_blank(),
    legend.text = element_text(face = 'bold', size = 12),
    legend.title = element_text(face = 'bold', size = 14),
    legend.position = c(0.8, 0.15),
  )

ggsave('plots/fig2_ideas/rel_abundance_core_families.png', 
       width = 8, height = 8)

core_family_rel_abund_table %>% 
  group_by(Family) %>% 
  summarize(mean = mean(rel_abund),
            sd = sd(rel_abund)) %>% 
  arrange(-mean)

  
  