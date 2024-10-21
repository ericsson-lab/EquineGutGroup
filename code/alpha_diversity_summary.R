library(rstatix)
library(randomForest)
library(pairwiseAdonis)
source('code/load_data.R')
source('code/data_beta_diversity.R')

main_color = '#416d3e'
  

## Figure 2A READ COUNTS
read_counts %>% 
  filter(!sampleid %in% c(samples_to_drop$sampleid, unknown_samples)) %>% 
  get_summary_stats(type = 'common') %>% 
  write_csv('stats/summary_stats/trimmed_read_counts.csv')

read_counts %>% 
  filter(!sampleid %in% c(samples_to_drop$sampleid, unknown_samples)) %>% 
  ggplot(aes(x = log10(trimmed_forward))) +
  geom_histogram(fill = main_color, color = main_color) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = scales::comma, breaks = c(0, 500, 1000, 1500)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, 6)) +
  labs(y = '# of Samples',
       x = "Log<sub>10</sub> Read Count") +
  theme(
    aspect.ratio = 2/1,
    axis.title.x = ggtext::element_markdown()
  )
ggsave('plots/sample_summary/forward_read_counts.png', width = 3.25,
       height = 3)

## Figure 2B FEATURE COUNTS
features_per_sample <- table %>% 
  pivot_longer(-featureid, names_to = 'sampleid', values_to = 'count') %>% 
  group_by(sampleid) %>% 
  summarize(count = sum(count)) %>% 
  filter(!sampleid %in% samples_to_drop$sampleid)

features_per_sample %>% 
  filter(sampleid %in% c(samples_to_drop$sampleid, unknown_samples))
# All samples removed

features_per_sample %>% 
  get_summary_stats(type = 'common') %>% 
  write_csv('stats/summary_stats/feature_counts.csv')

features_per_sample %>% 
  ggplot(aes(x = log10(count))) +
  geom_histogram(fill = main_color, color = main_color, bins = 30) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = scales::comma, breaks = c(0, 500, 1000, 1500)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)), 
                     limits = c(0, 6)) +
  labs(y = '# of\ Features',
       x = "Log<sub>10</sub> Feature Count") +
  theme(
    aspect.ratio = 2/1,
    axis.title.x = ggtext::element_markdown()
  )

ggsave('plots/sample_summary/feature_counts.png', width = 3.25,
       height = 3)


# Read in alpha diversity metrics
alpha_stats <- read_rds('data/240711_merged/data/alpha.stats.Rds') %>% 
  left_join(., metadata) %>% 
  rename(age = 'age_of_collection',
         # diagnosis = 'disease_category', 
         month = 'month_norm') %>% 
  mutate(age = as.double(age/365)) 

# Figure 2C Chao1 overall
chao1_summary <- alpha_stats %>% 
  summarise(mean = mean(chao1),
            sd = sd(chao1))

alpha_stats %>% 
  get_summary_stats(chao1, type = 'common') %>% 
  write_csv('stats/summary_stats/chao1_all_samples.csv')

alpha_stats %>% 
  arrange(chao1) %>% 
  mutate(n=row_number()) %>% 
  ggplot(aes(x = n, y = chao1)) +
  geom_rect(aes(ymin = chao1_summary$mean - chao1_summary$sd, 
                ymax = chao1_summary$mean + chao1_summary$sd, 
                xmin = min(after_stat(x)),
                xmax = max(after_stat(x))),
            fill = '#96b885', color = NA) +
  geom_point() +
  geom_hline(aes(yintercept = chao1_summary$mean),
             linetype = 2, linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     limits = c(0, 3000)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03)),
                     breaks = c(0, 800, 1600, 2400), limits = c(0, 2363)) +
  xlab('# of Samples') +
  ylab('Chao1 Index') +
  theme(
    aspect.ratio = 3/2,
    # axis.text.x = element_blank(),
    # axis.ticks.x= element_blank()
  )

ggsave('plots/alpha_stats/chao1/all_samples.png', width = 4, height = 3)

# Figure 2D Shannon data overall
shannon_summary <- alpha_stats %>% 
  summarise(mean = mean(shannon),
            sd = sd(shannon))
alpha_stats %>% 
  get_summary_stats(shannon, type = 'common') %>% 
  write_csv('stats/summary_stats/shannon_all_samples.csv')

alpha_stats %>% 
  arrange(shannon) %>% 
  mutate(n=row_number()) %>% 
  ggplot(aes(x = n, y = shannon)) +
  geom_rect(aes(ymin = shannon_summary$mean - shannon_summary$sd,
                ymax = shannon_summary$mean + shannon_summary$sd,
                xmin = min(after_stat(x)),
                xmax = max(after_stat(x))),
  fill = '#96b885', color = NA) +
  geom_point() +
  geom_hline(aes(yintercept = shannon_summary$mean),
             linetype = 2, linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     limits = c(0, 8)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03)),
                     breaks = c(0, 800, 1600, 2400), limits = c(0, 2400)) +
  xlab('# of Samples') +
  ylab('Shannon Index') +
  theme(
    aspect.ratio = 3/2,
    # axis.text.x = element_blank(),
    # axis.ticks.x= element_blank()
  )

ggsave('plots/alpha_stats/shannon/all_samples.png', width = 3, height = 3)

# Figure 2E
## Random forest
chao1_data <- alpha_stats %>% 
  select(sampleid, chao1, sex, breed_type, age, 
         country, diagnosis,month) %>% 
  drop_na() %>% 
  column_to_rownames(var = 'sampleid')
# 
# index = sample(2, nrow(chao1_data), replace = T, prob = c(0.8, 0.2))
# train = chao1_data[index == 1, ]
# test = chao1_data[index == 2, ]

chao1_rf_model <- randomForest(chao1 ~ ., data = chao1_data, proximity = T, 
                               importance = T, ntree = 1000)
print(chao1_rf_model)
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
  geom_col(fill = '#416d3e')+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, 150),
                     breaks = c(0, 75, 150)) +
  ggprism::theme_prism() +
  ylab('%IncMSE') +
  theme(
    aspect.ratio = 3/1.25,
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  coord_flip() 
ggsave('plots/alpha_stats/rf_importance_chao1.png', width = 3, height = 3)

shannon_data <- alpha_stats %>% 
  select(sampleid, shannon, sex, breed_type, age, month, 
         country, diagnosis) %>% 
  drop_na() %>%   
  column_to_rownames(var = 'sampleid')

shannon_rf_model <- randomForest(shannon ~ ., data = shannon_data, proximity = T, 
                                 importance = T, ntree = 1000)
plot(shannon_rf_model)
importance(shannon_rf_model)
print(shannon_rf_model)

shannon_importance <- importance(shannon_rf_model) %>% 
  as_tibble(rownames = 'category') %>% 
  select(category, `%IncMSE`) %>% 
  rename(IncMSE_shannon = `%IncMSE`) 

shannon_importance %>% 
  mutate(category = str_replace_all(category, '_', ' '),
         category = str_to_title(category)) %>% 
  ggplot(aes(x = factor(category, levels = c('Sex', 'Breed Type', 'Month', 
                                             'Country', 'Diagnosis', 'Age')), 
             y = IncMSE_shannon)) +
  geom_col(fill = '#416d3e') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, 150),
                     breaks = c(0, 75, 150)) +
  ggprism::theme_prism() +
  ylab('%IncMSE') +
  theme(
    aspect.ratio = 3/1.25,
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  coord_flip() 

ggsave('plots/alpha_stats/rf_importance_shannon.png', width = 3, height = 3)


chao1_importance %>% 
  left_join(., shannon_importance) %>% 
  write_csv('stats/summary_stats/alpha_stat_rd_incMSE.csv')


# Supplementary Figure 3A-B
plot_alpha_age <- function(metric, ylab, upper_limit){
  alpha_stats %>% 
    select(chao1, shannon, age) %>% 
    drop_na() %>% 
    ggplot(aes(x = age, y = {{metric}}, color = age)) +
    geom_point() +
    scale_color_gradient(low = '#ecf0ec', high = '#416d3e',
                         limits = c(0, 40), name = 'Age') +
    scale_y_continuous(limits = c(0, upper_limit), expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(limits = c(0, 40), expand = expansion(mult = c(0, 0.05))) +
    
    labs(y = ylab,
         x = 'Age (Years)') +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      aspect.ratio = 1/1.75,
      legend.position = 'none'
    )
}

alpha_div_data <- alpha_stats %>% 
  select(chao1, shannon, age) %>% 
  drop_na()

# Fail noramility, using spearman correlation
alpha_div_data %>% 
  shapiro_test(chao1, shannon)

# chao1 (supplementary figure 3a)
plot_alpha_age(chao1, 'Chao1 Index', 3000)
ggsave('plots/alpha_stats/chao1_age.png', width = 6, height = 4)
alpha_div_data %>% 
  cor_test(age, chao1, method = 'spearman')
# shannon (supplementary figure 3b)
plot_alpha_age(shannon, 'Shannon Index', 7) +
  scale_y_continuous(breaks = seq(0,7,1), expand = expansion(mult = c(0, 0.1)))
alpha_div_data %>% 
  cor_test(age, shannon, method = 'spearman')
ggsave('plots/alpha_stats/shannon_age.png', width = 6, height = 4)



# tictoc::tic()
# bc_dist <- generate_dist(table = table, distance = 'bray')
# bc_dist %>%
#   write_rds('data/240711_merged/data/bc_dist_raw.Rds')
# tictoc::toc()
pcoa_metadata <- chao1_data %>% 
  rownames_to_column(var = 'sampleid') %>% 
  as_tibble()

bc_dist <- read_rds('data/240711_merged/data/bc_dist_raw.Rds') %>% 
  usedist::dist_subset(., idx = pcoa_metadata$sampleid)

bc_pcoa <- generate_pcoa(dist = bc_dist, metadata_df = pcoa_metadata)

pcoa_plot <- bc_pcoa[[1]]
pcoa_plot_axis <- bc_pcoa[[2]]

# tictoc::tic()
# overall_adonis_res <- adonis2(bc_dist ~ age + diagnosis + country + breed_type + month + sex,
#                permutations = 9999, data = pcoa_metadata)
# overall_adonis_res %>% 
#   as.tibble(rownames = 'effect') %>% 
#   write_csv('stats/summary_stats/bc_dist_anova_res_all.csv')
# tictoc::toc()

overall_adonis_res %>% 
  as.tibble(rownames = 'effect') %>% 
  drop_na() %>% 
  mutate(effect = factor(effect, levels = c('sex',  'breed_type', 'month',
                                            'country', 'diagnosis', 'age'))) %>% 
  ggplot(aes(x = effect, y = R2))+
  geom_col(fill = '#416d3e') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, 0.05),
                     breaks = c(0, 0.05), 
                     labels = scales::percent
                     ) +
  ggprism::theme_prism() +
  ylab('R^2^') +
  theme(
    aspect.ratio = 3/1.25,
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = ggtext::element_markdown()
  ) +
  coord_flip() 
 ggsave('plots/alpha_stats/anova_res_effect_size.png', width = 3, height = 3)

 # Figure 2f
pcoa_plot %>% 
  ggplot(aes(x = PCo1, y = PCo2, color = age)) +
  geom_point(size = 3, alpha = 1) +
  scale_color_gradient(low = '#ecf0ec', high = '#416d3e',
                       limits = c(0, 40), name = 'Age\n(years)') +
  xlab(glue::glue('PCo1 - {round(pcoa_plot_axis[1], 2)}%')) +
  ylab(glue::glue('PCo2 - {round(pcoa_plot_axis[2], 2)}%')) +
  theme(
    aspect.ratio = 1/1.25,
    legend.ticks = element_blank(),
    legend.text = element_text(face = 'bold', size = 12),
    legend.title = element_text(face = 'bold', size = 14, hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  )
ggsave('plots/sample_summary/overall_bc_pcoa_age.png', width = 6, height = 4)


# Supplementary Figure 3E
diagnosis_pcoa_plot <- pcoa_plot %>%
  mutate(diagnosis = fct_lump(diagnosis, 3, other_level = 'Other')) 

diagnosis_pcoa_plot %>% 
  ggplot(aes(x = PCo1, y = PCo2, 
             color = factor(diagnosis, levels = c('Healthy', 'Gastrointestinal (Acute)',
                                                  'Musculoskeletal', 'Other')))) +
  geom_point(size = 3, alpha = 0.6,
             data = diagnosis_pcoa_plot %>% filter(diagnosis == 'Other'))+
  geom_point(size = 3, alpha = 0.6,
             data = diagnosis_pcoa_plot %>% filter(diagnosis != 'Other')) +
  scale_color_manual(values = c('blue', 'red', main_color, 'gray'),
                     labels = c('Healthy', 'Acute GI', 'Musculoskeletal', 'Other'),
                     breaks = c('Healthy', 'Gastrointestinal (Acute)', 'Musculoskeletal', 'Other'),
                     name = 'Disease\nCategory') +
  xlab(glue::glue('PCo1 - {round(pcoa_plot_axis[1], 2)}%')) +
  ylab(glue::glue('PCo2 - {round(pcoa_plot_axis[2], 2)}%')) +
  theme(
    aspect.ratio = 1/1,
    legend.ticks = element_blank(),
    legend.text = element_text(face = 'bold', size = 12),
    legend.title = element_text(face = 'bold', size = 14, hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1))
  )
ggsave('plots/sample_summary/overall_bc_pcoa_diagnosis.png', width = 6, height = 4)

pairwise.adonis2(bc_dist ~ country, data =diagnosis_pcoa_plot, nperm = 9999)


# Supplementary Figure 3h
pcoa_plot %>% 
  ggplot(aes(x = PCo1, y = PCo2, 
             color = factor(country, levels = c('USA', 'UK', 'ASTL')))) +
  geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = c('blue','red', main_color), 
                     labels = c('USA', 'UK', 'AUS'),
                     name = 'Country') +
  xlab(glue::glue('PCo1 - {round(pcoa_plot_axis[1], 2)}%')) +
  ylab(glue::glue('PCo2 - {round(pcoa_plot_axis[2], 2)}%')) +
  theme(
    aspect.ratio = 1/1,
    legend.ticks = element_blank(),
    legend.text = element_text(face = 'bold', size = 12),
    legend.title = element_text(face = 'bold', size = 14, hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1))
  )
ggsave('plots/sample_summary/overall_bc_pcoa_country.png', width = 6, height = 4)

# Supplementary Figure 3f
alpha_stats %>% 
  mutate(country = factor(country, levels = c('USA', 'UK', 'ASTL'))) %>% 
  ggplot(aes(x = fct_rev(country), y = chao1, color = country)) +
  geom_point(position = position_jitter()) +
  stat_summary(geom = 'bar', fun = 'mean', fill = NA, 
               linewidth = 1, color = 'black') +
  stat_summary(geom = 'errorbar',
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               width = 0.5, color = 'black', 
               linewidth = 1) +
  scale_color_manual(values = c('blue', 'red', main_color),
                     labels = c('USA', 'UK', 'AUS')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     breaks = c(0, 1500, 3000)) +
  ylab('Chao1 Index') +
  theme(
    aspect.ratio = 3/2,
    legend.position = 'none',
    axis.title.y = element_blank(),
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
    ) +
  coord_flip()

ggsave('plots/sample_summary/chao1_country.png', width = 5, height = 4)

# Supplementary Figure 3g
alpha_stats %>% 
  mutate(country = factor(country, levels = c('USA', 'UK', 'ASTL'))) %>% 
  ggplot(aes(x = fct_rev(country), y = shannon, color = country)) +
  geom_point(position = position_jitter()) +
  stat_summary(geom = 'bar', fun = 'mean', fill = NA, 
               linewidth = 1, color = 'black', ) +
  stat_summary(geom = 'errorbar',
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               width = 0.5, color = 'black', 
               linewidth = 1) +
  scale_color_manual(values = c('blue', 'red', main_color),
                     labels = c('USA', 'UK', 'AUS'),
                     breaks = c('USA', 'UK', 'ASTL')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.23)),
                     breaks = seq(0,8,2)) +
  ylab('Shannon Index') +
  theme(
    aspect.ratio = 3/2,
    legend.position = 'none',
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 16)
  ) +
  coord_flip()
ggsave('plots/sample_summary/shannon_country.png', width = 5, height = 4)



alpha_stats %>% 
  kruskal_test(chao1 ~ country)
alpha_stats %>% 
  pairwise_wilcox_test(chao1 ~ country, p.adjust.method = 'BH')
alpha_stats %>% 
  kruskal_test(shannon ~ country)
alpha_stats %>% 
  pairwise_wilcox_test(shannon ~ country, p.adjust.method = 'BH')



alpha_stats %>% 
  mutate(diagnosis = case_when(diagnosis == 'Gastrointestinal (Acute)' ~ 'Acute GI',
                               TRUE ~ diagnosis)) %>% 
  ggplot(aes(y = chao1, 
             x = factor(fct_lump(diagnosis, 3), 
                        levels = rev(c('Healthy', 'Acute GI',
                                       'Musculoskeletal', 'Other'))),  
             color = factor(fct_lump(diagnosis, 3), 
                            levels = rev(c('Healthy', 'Acute GI',
                                           'Musculoskeletal', 'Other'))))) +
  geom_point(position = position_jitter()) +
  stat_summary(geom = 'bar', fun = 'mean', fill = NA, 
               linewidth = 1, color = 'black') +
  stat_summary(geom = 'errorbar',
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               width = 0.5, color = 'black', 
               linewidth = 1)  +
  scale_color_manual(values = rev(c('blue', 'red', main_color, 'lightgray')),
                     breaks = rev(c('Healthy', 'Acute GI', 
                                    'Musculoskeletal', 'Other')),
                     labels = c('Healthy', 'Acute GI', 'Musculoskeletal', 'Other')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.23)),
                     breaks = c(0, 1500, 3000)) +
  ylab('Chao1 Index') +
  theme(
    aspect.ratio = 3/2,
    legend.position = 'none',
    axis.title.y = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16)
  ) +
  coord_flip()
ggsave('plots/sample_summary/chao1_disease_category.png', width = 5, height = 4)

alpha_stats %>% 
  mutate(diagnosis = case_when(diagnosis == 'Gastrointestinal (Acute)' ~ 'Acute GI',
                               TRUE ~ diagnosis)) %>% 
  ggplot(aes(y = shannon, 
             x = factor(fct_lump(diagnosis, 3), 
                        levels = rev(c('Healthy', 'Acute GI',
                                       'Musculoskeletal', 'Other'))),  
             color = factor(fct_lump(diagnosis, 3), 
                            levels = rev(c('Healthy', 'Acute GI',
                                           'Musculoskeletal', 'Other'))))) +
  geom_point(position = position_jitter()) +
  stat_summary(geom = 'bar', fun = 'mean', fill = NA, 
               linewidth = 1, color = 'black') +
  stat_summary(geom = 'errorbar',
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               width = 0.5, color = 'black', 
               linewidth = 1) +
  scale_color_manual(values = rev(c('blue', 'red', main_color, 'lightgray')),
                     breaks = rev(c('Healthy', 'Acute GI', 
                                    'Musculoskeletal', 'Other')),
                     labels = c('Healthy', 'Acute GI', 'Musculoskeletal', 'Other')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.23)),
                     breaks = seq(0,8,2)) +
  ylab('Chao1 Index') +
  theme(
    aspect.ratio = 3/2,
    legend.position = 'none',
    axis.title.y = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16)
  ) +
  coord_flip()


alpha_stats %>% 
  kruskal_test(chao1 ~ diagnosis)
alpha_stats %>% 
  pairwise_wilcox_test(chao1 ~ diagnosis, p.adjust.method = 'BH') %>% 
  write_csv('stats/alpha_stats/chao1_primarydiagnosis_wilcoxon_res.csv')
alpha_stats %>% 
  kruskal_test(shannon ~ diagnosis)
alpha_stats %>% 
  pairwise_wilcox_test(shannon ~ diagnosis, p.adjust.method = 'BH') %>% 
  write_csv('stats/alpha_stats/shannon_primarydiagnosis_wilcoxon_res.csv')


pairwise.adonis2(bc_dist ~ diagnosis, data = pcoa_plot, nperm = 9999 )


university_pal <- c(
  'Auburn University' = "#9E0142",
  'Louisiana State University' = "#D53E4F",
  'North Carolina State University' = "#F46D43",
  'University of Minnesota' = "#FDAE61",
  'University of Missouri' = "#ABDDA4" ,
  'University of Liverpool' = "#66C2A5",
  'University of Queensland' = "#3288BD",
  'University of Sydney' = "#5E4FA2"
)
alpha_stats %>% 
  ggplot(aes(x = fct_relevel(university, 'University of Liverpool',
                             after = 5), 
             y = chao1, color = university)) +
  geom_point(position = position_jitter(), size = 2) +
  stat_summary(geom = 'bar', fun = 'mean',
               color = 'black', linewidth = 1,
               fill = NA) +
  stat_summary(geom = 'errorbar', color = 'black', width = 0.5,
               linewidth = 1,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_color_manual(values = university_pal) +
  coord_flip() +
  ylab('Chao1 Index') +
  theme(
    aspect.ratio = 2/1,
    axis.title.y = element_blank(),
    legend.position = 'none'
  )

ggsave('plots/alpha_stats/chao1_university.png', width = 6, height = 5)

alpha_stats %>% 
  kruskal_test(chao1 ~ university)
alpha_stats %>% 
  pairwise_wilcox_test(chao1 ~ university, 
                       p.adjust.method = 'BH') %>% 
  write_csv('stats/summary_stats/chao1_university_wilcoxon.csv')
alpha_stats %>% 
  pairwise_wilcox_test(shannon ~ university, 
                       p.adjust.method = 'BH') %>% 
  write_csv('stats/summary_stats/shannon_university_wilcoxon.csv')

library(pairwiseAdonis)
pairwise.adonis2(bc_dist ~ university, 
                 data = alpha_stats %>% , nperm = 9999)


alpha_stats %>% 
  ggplot(aes(x = fct_relevel(university, 'University of Liverpool',
                             after = 5), 
             y = shannon, color = university)) +
  geom_point(position = position_jitter(), size = 2) +
  stat_summary(geom = 'bar', fun = 'mean',
               color = 'black', linewidth = 1,
               fill = NA) +
  stat_summary(geom = 'errorbar', color = 'black', width = 0.5,
               linewidth = 1,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_color_manual(values = university_pal) +
  coord_flip() +
  ylab('Shannon Index') +
  theme(
    aspect.ratio = 2/1,
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )
ggsave('plots/alpha_stats/shannon_university.png', width = 6, height = 5)


pcoa_plot <- pcoa_plot %>% 
  left_join(., metadata) 
pcoa_plot %>% 
  ggplot(aes(x = PCo1, y = PCo2, color = university)) +
  geom_point(size = 2,
             data = pcoa_plot %>% 
               filter(university == 'University of Missouri')) +
  geom_point(size = 2,
             data = pcoa_plot %>% 
               filter(university != 'University of Missouri')) +
  scale_color_manual(values = university_pal, name = 'University') +
  labs(
    x = glue::glue('PCo1 - {round(pcoa_plot_axis[1], 2)}%'),
    y = glue::glue('PCo2 - {round(pcoa_plot_axis[2], 2)}%')
  ) +
  theme(
    aspect.ratio = 2/1.5,
    legend.text = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold', size = 14, hjust = 0.5),
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4))
  )
ggsave('plots/alpha_stats/pcoa_university.png', width = 12, height = 5)


bc_dist <- read_rds('data/240711_merged/data/bc_dist_raw.Rds') %>% 
  usedist::dist_subset(., idx = alpha_stats$sampleid)
pairwise.adonis2(bc_dist ~ university, nperm = 9999, data = alpha_stats)



