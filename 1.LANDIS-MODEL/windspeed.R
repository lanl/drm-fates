library(here)
library(tidyverse)
library(ggthemes)

ws <- read_csv(here("1.LANDIS-MODEL","QF_outputs","windspeeds.csv"))

ws_long <- ws %>%
  pivot_longer(cols = 3:6,
               names_to = "var",
               values_to = "val")

fire_colors <- solarized_pal()(7)[c(7,5,6,3)]
ws_long %>% 
  filter(var == "midcanopy_windspeed") %>%
  filter(timestep < 70) %>%
  mutate(cycle = factor(cycle, labels = c("Fire 1","Fire 2","Fire 3","Fire 4"))) %>%
  ggplot() +
  geom_line(aes(timestep*30, val, color = cycle)) +
  scale_color_manual(values = fire_colors) +
  labs(x="Time (s)",
       y="Subcanopy wind speed (m/s)") +
  theme_bw() +
    theme(legend.title = element_blank())
