library(tidyverse)
library(janitor)
library(readxl)
library(maps)
library(usdata)
library(ggstar)

source('code/load_data.R')

usa <- map_data("state")

# 3 horses did not have location data in US
horse_counts_by_state <- metadata %>% 
  filter(country == 'USA') %>% 
  select(horse, state) %>% 
  drop_na() %>% 
  distinct()

state_counts <- horse_counts_by_state %>% 
  count(state) %>% 
  mutate(state = abbr2state(state),
         state = tolower(state)) %>% 
  rename(region = 'state')


schools_lat_long <- tibble(
  country = c('US', 'US', 'US', 'US', 'US', 'UK', 'ASTL', 'ASTL'),
  school = c('University of Missouri', 'University of Minnesota', 'Lousiana State University', 'Aurburn University', 'North Carolina State University', 'University of Liverpool', 'University of Sydney', 'University of Queensland'),
  x = c( -92.31806, -93.18338, -91.19343, -85.51144,-78.70435, -3.02733, 150.65821, 152.33618),
  y = c( 38.94007, 44.98183, 30.41369, 32.58954, 35.79911,  53.28854, -34.03039, -27.55052)
) 

midwest = c('north dakota', 'south dakota', 'nebraska', 'kansas', 'missouri', 'iowa', 'minnesota', 'wisconsin', 'illinois', 'indiana', 'michigan', 'ohio')
southeast = c('arkansas', 'louisiana', 'mississippi', 'alabama', 'georgia', 'florida', 'south carolina', 'north carolina', 'virginia', 'west virginia', 'kentucky', 'tennessee')
usa %>%
  filter(region %in% c(midwest, southeast, 'oklahoma')) %>% 
  left_join(., state_counts) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(aes(group = group, fill = n), 
               color = 'white') +
  geom_star(aes(x = x, y = y), inherit.aes = F, data = schools_lat_long %>% filter(country == 'US'),
            fill = 'white', color = 'black', size = 3, starstroke = 1) +
  scale_fill_gradient(low = '#96b885', high = '#3a502f',
                      na.value = 'gray', limits = c(0, 750),
                      breaks = c(1, 250, 500, 750),
                      name = 'Number\nof Patients') +
  theme_void() +
  coord_map() +
  theme(legend.title = element_text(face = 'bold', hjust = 0.5),
        legend.text = element_text(face = 'bold', size = 12),
        legend.direction =  'horizontal', 
        legend.position = c(0.9, 0.8),
        legend.position = 'none'
        ) +
  guides(
    fill = guide_colorbar(ticks.colour = NA, title.position = 'top')
  )
ggsave('plots/maps/midwest_southeast_ok.png', width = 8, height = 4)


australia_map <- map_data(map = 'world', region = "Australia")
uk_map <- map_data(map = 'world', region = "UK")

# count patients in UK and ASTL
metadata %>% 
  filter(country != 'USA') %>% 
  select(horse, country) %>% 
  distinct() %>%
  count(country)
# ASTL 137
# UK   178

uk_map %>% 
  mutate(fill = 178) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(aes(group = group, color = subregion, fill = fill), color = 'white')+
  geom_star(aes(x = x, y = y), inherit.aes = F, data = schools_lat_long %>% filter(country == 'UK'),
            fill = 'white', color = 'black', size = 6, starstroke = 2) +
  scale_fill_gradient(low = '#96b885', high = '#3a502f',
                      na.value = 'gray', limits = c(0, 750)) +
  theme_void() +
  coord_map() +
  theme(
    legend.position = 'none'
  )
ggsave('plots/maps/uk_map.png', width = 4, height = 4)


australia_map %>% 
  mutate(fill = 137) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(aes(group = group, fill = fill), color = 'white')+
  geom_star(aes(x = x, y = y), inherit.aes = F, data = schools_lat_long %>% filter(country == 'ASTL'),
            fill = 'white', color = 'black', size = 4, starstroke = 1.5) +
  scale_fill_gradient(low = '#96b885', high = '#3a502f',
                      na.value = 'gray', limits = c(0, 750)) +
  theme_void() +
  coord_map()  +
  theme(legend.position = 'none')

ggsave('plots/maps/astl_map.png', width = 4, height = 4)



## Supplementary Figure 1
mo_county_map_data <- map_data("county") %>% 
  filter(region == 'missouri')

mo_horse_zip <- metadata %>% 
  filter(state == 'MO') %>% 
  select(horse, state, zip) %>% 
  distinct() %>% 
  drop_na() %>% 
  rename(zip_code = 'zip') %>% 
  mutate(zip_code = as.double(zip_code))

mo_zip_names <- read_csv('data/Missouri_Zip_Codes_by_County_City.csv') %>% 
  clean_names()

setdiff(mo_horse_zip$zip_code, mo_zip_names$zip_code) 
# all zip codes from horse data are in mo data

mo_detected_counties <- mo_zip_names %>% 
  filter(zip_code %in% mo_horse_zip$zip_code) %>% 
  select(county, zip_code) %>% 
  distinct()

mo_county_counts <- mo_horse_zip %>% 
  left_join(., mo_detected_counties) %>% 
  group_by(county) %>% 
  count() %>% 
  mutate(county = tolower(county)) %>% 
  rename(subregion = 'county')

counties_covered <- mo_county_counts$subregion %>% unique() %>% length()
mo_counties <- mo_zip_names$county %>% unique() %>% length()
(counties_covered / mo_counties) * 100

mo_county_map_data %>% 
  left_join(., mo_county_counts) %>% 
  ggplot(aes(y = lat, x = long)) +
  geom_polygon(aes(group = subregion, fill = n), 
               color = 'white') +
  geom_star(aes(x = x, y = y), inherit.aes = F, 
            data = schools_lat_long %>% filter(school == 'University of Missouri'), color = 'black', fill = 'white', 
            size = 3) +
  coord_map() +
  scale_fill_gradient(low = '#96b885', high = '#3a502f',
                      na.value = 'gray', limits = c(0, 150),
                      breaks = c(1, 50, 100, 150),
                      name = '# of\nPatients') +
  theme_void() +
  theme(legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', hjust = 0.5, size = 14),
        legend.direction = 'horizontal',
        legend.position = c(0.8, 0.85)) +
  guides(
    fill = guide_colorbar(ticks.colour = NA, title.position = 'top')
  )

ggsave('plots/maps/mo_map.png', width = 6, height = 6, bg = 'white')




