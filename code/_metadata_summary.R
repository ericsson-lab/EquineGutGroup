source('code/load_data.R')
library(rstatix)

theme_set(ggprism::theme_prism(base_size = 18))

main_color = '#416d3e'
light_color =  '#92bf8f'
    
# Figure 1B, samples per year
# 49 samples do not have a collection date
metadata %>% 
  select(sampleid, collection_date) %>% 
  ggplot(aes(x = collection_date)) +
  geom_histogram(fill = '#416d3e', color  = '#416d3e') +
  scale_y_continuous(breaks = c(0, 200, 400), expand = expansion(mult = c(0, 0.05))) +
  ylab('# of Samples') +
  theme(
    aspect.ratio = 3/1,
    legend.title = element_text(face = 'bold'),
    legend.text = element_text(face = 'bold'),
    legend.position = c(1.1, 0.9),
    axis.title.y = element_blank()
  ) +
  coord_flip()
ggsave('plots/sample_summary/collection_year.png', height = 5, width = 6)

metadata %>% 
  drop_na(collection_date) %>% 
  summarize(min = min(collection_date),
            max = max(collection_date))

# Figure 1C, samples per year
# 49 samples do not have a collection date
metadata %>% 
  select(country, collection_date) %>%
  mutate(month = month(collection_date),
         hemisphere = case_when(country == 'ASTL' ~ 'S',
                                TRUE ~ 'N')) %>% 
  group_by(month, hemisphere) %>% 
  drop_na() %>%
  count() %>% 
  ggplot(aes(x= fct_rev(as.factor(month)), y = n, fill = hemisphere))  +
  geom_col() +
  scale_y_continuous(breaks = c(0, 200, 400), expand = expansion(mult = c(0, 0.05)))  +
  scale_x_discrete(labels = c(rev(month.name))) +
  scale_fill_manual(values = c('#416d3e', '#92bf8f'),
                    name = 'Hemi') +
  ylab('# of Samples') +
  theme(
    aspect.ratio = 3/1,
    legend.title = element_text(face = 'bold', 
                                angle = 90,
                                hjust = 0.5),
    legend.text = element_text(face = 'bold', size = 16),
    legend.box = 'vertical',
    legend.title.position = 'left',
    legend.position = c(0.8, 0.9),
    axis.title.y = element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 15))
  ) +
  coord_flip()

ggsave('plots/sample_summary/time_of_year.png', height = 5, width = 6)


# Figure 1D, equine_sex
metadata %>%
  mutate(equine_sex = str_to_sentence(equine_sex)) %>% 
  ungroup() %>% 
  select(horse, equine_sex, sex) %>% 
  distinct() %>%
  filter(sex != 'Unknown') %>% 
  group_by(equine_sex, sex) %>% 
  count() %>% 
  ggplot(aes(x = fct_reorder(equine_sex, n), y = n)) +
  geom_col(fill = '#416d3e') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = c(0, 300, 600),
                     limits = c(0, 650),
                     name = '# of Horses') +
  facet_wrap(~sex, scales = 'free', nrow = 2, strip.position = 'left') +
  coord_flip() +
  theme(
    aspect.ratio = 3/2, 
    axis.title.y = element_blank(),
    strip.text = element_text(size = 20),
    strip.placement = 'outside'
  )

ggsave('plots/sample_summary/sex_distribution.png', width = 5, height = 6)

metadata %>%
  mutate(equine_sex = str_to_sentence(equine_sex)) %>% 
  ungroup() %>% 
  select(horse, equine_sex, sex) %>% 
  distinct() %>%
  filter(sex != 'Unknown') %>%
  group_by(equine_sex, sex) %>% 
  count() %>% 
  arrange(sex, -n) %>% 
  write_csv('data/240711_merged/summary/equine_sex.csv')

metadata %>%
  mutate(equine_sex = str_to_sentence(equine_sex)) %>% 
  ungroup() %>% 
  select(horse, sex) %>% 
  distinct() %>%
  count(sex) %>% 
  mutate(p = n / sum(n))


# Figure 1E Age at collection
# 60 with no age
metadata %>% 
  mutate(age_of_collection = as.double(age_of_collection)/365) %>% 
  ggplot(aes(x = age_of_collection)) +
  geom_histogram(color = '#416d3e', fill = '#416d3e', bins = 30, width = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks = c(0, 100, 200)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = '# of Samples',
       x = 'Sampling Age (years)') +
  theme(
    aspect.ratio = 3/1
  ) +
  coord_flip()
ggsave('plots/sample_summary/age_of_collection.png', height = 6, width = 5)

metadata %>% 
  mutate(age_of_collection = as.double(age_of_collection)/365) %>% 
  drop_na(age_of_collection) %>% 
  get_summary_stats(age_of_collection, type = 'common')
  
# Figure 1F Breeds
metadata %>% 
  mutate(breed = str_to_title(breed)) %>% 
  ungroup() %>% 
  select(horse, breed) %>% 
  distinct() %>% 
  count(breed) %>% 
  arrange(-n) %>% 
  mutate(breed = case_when(breed == 'Missouri Fox Trotter' ~ 'MO Fox Trotter',
                           TRUE ~ breed),
         row_number = row_number(),
         other = row_number > 11,
         group = ifelse(other, 'Other', breed),
         group = case_when(group == 'Unknown' ~ 'Other',
                           TRUE ~ group)) %>% 
  summarize(.by = 'group', n = sum(n)) %>% 
  ggplot(aes(x = fct_relevel(fct_reorder(group, n), 'Other', after = 0), y = n)) + 
  geom_col(fill = '#416d3e') +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)),
                     breaks = c(0, 200, 400),
                     limits = c(0, 400)) +
  coord_flip() +
  ylab('# of Horses') +
  theme(
    aspect.ratio = 3/1,
    axis.title.y = element_blank()
  )

ggsave('plots/sample_summary/breed_counts.png',  width = 5, height = 6)  

metadata %>% 
  mutate(breed = str_to_title(breed)) %>% 
  ungroup() %>% 
  select(horse, breed) %>% 
  distinct() %>% 
  count(breed) %>% 
  arrange(-n) %>% 
  print(n = 100) %>% 
  write_csv('data/240711_merged/summary/breed_counts_horses.csv')


# Supplementary Figure 2a 
# Breed_count
metadata %>% 
  select(horse, breed_type) %>% 
  distinct() %>% 
  count(breed_type) %>% 
  arrange(-n) %>% 
  ggplot(aes(x = fct_relevel(fct_reorder(breed_type, n), 'Other', after = 0), y = n)) + 
  geom_col(fill = main_color) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)),
                     breaks = c(0, 200, 400),
                     limits = c(0, 400)) +
  coord_flip() +
  ylab('# of Patients') +
  theme(
    aspect.ratio = 3/1,
    axis.title.y = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  )
ggsave(
  'plots/sample_summary/breed_type.png', width = 6, height = 5
)

metadata %>% 
  select(horse, breed_type) %>% 
  distinct() %>% 
  count(breed_type) %>% 
  arrange(-n)

# Figure 1G Sample & Host Counts
## DONE
horses_count <- metadata %>% 
  ungroup() %>% 
  select(horse, visit, disease_category) %>% 
  distinct() %>% 
  count(disease_category) %>% 
  rename(n_horses = 'n')

horses_count %>% 
  write_csv('data/240711_merged/summary/horses_disease_category.csv')
sample_count <- metadata %>% 
  ungroup() %>% 
  count(disease_category)
sample_count %>% 
  write_csv('data/240711_merged/summary/sample_disease_category.csv')


sample_count %>% 
  left_join(., horses_count) %>% 
  drop_na() %>% 
  filter(disease_category != 'Unknown') %>% 
  pivot_longer(-disease_category, values_to = 'n', names_to = 'group') %>% 
  ggplot(aes(x = fct_reorder(disease_category, n), y = n,
             fill = group), color = NA) +
  geom_col(position = position_dodge(-0.7), width = 0.7)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), breaks = c(0, 400, 800)) +
  scale_fill_manual(values = c('#416d3e', '#92bf8f'),
                    labels = c('Samples', 'Patients')) +
  labs(y = 'Count') +
  theme(
    axis.title.y = element_blank(),
    aspect.ratio = 3/1,
    legend.position = c(0.65, 0.1),
    legend.text = element_text(face = 'bold', size = 16)
  ) +
  coord_flip()

ggsave('plots/sample_summary/disease_category_sample_count.png', width = 5, height = 6)



## Supplementary Figure 2B
## All diagnoses
metadata %>% 
  select(horse, contains('disease_category')) %>% 
  distinct() %>% 
  pivot_longer(-c(horse, disease_category)) %>% 
  filter(value > 0) %>%
  ggplot(aes(x = fct_rev(fct_reorder2(horse, disease_category, fct_infreq(name))), 
             y = fct_rev(fct_infreq(name)), fill = value)) +
  geom_tile(width = 1, height = 1, show.legend = F,
            color = main_color) +
  scale_fill_manual(values = main_color) +
  scale_y_discrete(
    labels = rev(c('Gastrointestinal (Acute)', 'Musculoskeletal', 'Healthy',
                   'Ocular', 'Integumentary', 'Respiratory',
                   'Gastrointestinal (Chronic)','Neurologic',
                   'Reproductive', 'Cardiovascular', 'Dental',
                   'Systemic', 'Urinary','Immune-Mediated')),
  ) +
  xlab('Equine Patient') +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
  )

ggsave('plots/sample_summary/diagnosis.png', width = 7, height = 4)


# Supplementary Figure 2C
# Longitudinal data
metadata %>% 
  select(horse, visit, sample_order, disease_category) %>% 
  mutate(disease_category = case_when(disease_category == 'Healthy' ~ 'Healthy',
                                      TRUE ~ 'Clinical')) %>% 
  group_by(horse, visit, disease_category) %>% 
  count() %>% 
  group_by(disease_category, n) %>% 
  count(name = 'patients') %>% 
  ggplot(aes(x = n, y = (patients), color = fct_rev(disease_category))) +
  geom_line(linewidth = 2) + 
  scale_color_manual(values = c(main_color, light_color))+
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
  scale_x_continuous(breaks = seq(0,24,12)) +
  xlab('Samples / Patient') +
  ylab('# of Patients') +
  theme(
    aspect.ratio = 3/1,
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16,
                               face = 'bold'),
    legend.position = c(0.6, 0.8)) +
  guides(
    color = guide_legend(override.aes = list(shape = 14))
  )

ggsave(
  'plots/sample_summary/longitudinal_samples.png', width = 6, height = 5
)

## Euthanized patients 
metadata %>% 
  select(horse, euthanasia_death) %>% 
  mutate(euthanasia_death = str_to_lower(euthanasia_death)) %>% 
  distinct() %>% 
  count(euthanasia_death) %>% 
  mutate(p = n/sum(n))
