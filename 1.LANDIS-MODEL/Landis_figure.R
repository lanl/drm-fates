library(here)
library(tidyverse)
library(terra)
library(ggthemes)
library(tidyterra)
library(collapse)
library(patchwork)
library(sf)
library(gridExtra)
library(ggpattern)
library(scales)
library(ggh4x)

theme_new <- function(base_size = 9){
  theme_bw(base_size = base_size) %+replace%
    theme(
      axis.line.x = element_line(color="black", linewidth = 0.25),
      axis.line.y = element_line(color="black",linewidth = 0.25),
      axis.title = element_text(size = base_size),
      axis.text = element_text(colour="black", size=base_size-1),
      legend.key=element_rect(colour=NA, fill =NA),
      legend.text = element_text(colour="black", size=base_size-1),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(fill = NA, colour = NA),
      panel.background = element_rect(fill = "white", colour = "black"), 
      strip.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(size = base_size, margin = margin(4.4,4.4,4.4,4.4)),
    )
}
theme_set(theme_new())
col1 = 3.5
col2 = 7.16


read_perimeters <- function(cyc){
  perim <- rast(here("1.LANDIS-MODEL","QF_outputs",paste0("fuel_dens_final_cycle",cyc,".tif")))
  perim <- as.polygons(perim)
  names(perim) <- "IN"
  perim <- perim[perim$IN==1]
  perim$cycle <- cyc
  return(perim)
}

OC <- rast(here("1.LANDIS-MODEL","LANDIS_run","output-community-cycle0_cropped.tif"))
origin <- ext(OC)[c(1,3)]

FIA <- read_csv(here("9.FIA","FIA_raw","REF_SPECIES.csv"))

postlandis <- list()
for(i in 1:5){
  temp <- read_csv(here("1.LANDIS-MODEL","LANDIS_run",paste0("Treelist_postlandis_cycle",i-1,".csv")))
  temp$cycle <- i-1
  temp$stage <- "postlandis"
  postlandis[[i]] <- temp
}
postlandis <- bind_rows(postlandis)

postfire <- list()
for(i in 1:4){
  temp <- read_csv(here("1.LANDIS-MODEL","LANDIS_run",paste0("Treelist_postfire_cycle",i,".csv")))
  temp$cycle <- i
  temp$stage <- "postfire"
  postfire[[i]] <- temp
}
postfire <- bind_rows(postfire)

alldata <- bind_rows(postlandis,postfire)

OC_df <- as.data.frame(OC, xy=T)
names(OC_df) <- c("x","y","MapCode")

spgrp <- FIA %>%
  select(SPECIES_SYMBOL,COMMON_NAME,MAJOR_SPGRPCD,GENUS,SPECIES)

alldata_spgrp <- left_join(alldata, spgrp, by = join_by("SPECIES_SYMBOL")) %>%
  mutate(GENUS = if_else(GENUS=="Lithocarpus","Notholithocarpus",GENUS)) %>%
  mutate(sci_name = paste(GENUS,SPECIES),
         func_grp = if_else(MAJOR_SPGRPCD %in% c(1,2), "Conifer", "Hardwood"))

#create all possible bins
all_x_bins <- list()
i <- 1
for(bin in factor(levels(cut(0, breaks = seq(0,750,150))))){
  all_x_bins[[i]] <- rep(bin, length(factor(levels(cut(0, breaks = seq(0,750,150))))))
  i <- i+1
}
all_x_bins<-unlist(all_x_bins)
all_bins <- data.frame(X_bins = rep(factor(all_x_bins),5),
                       Y_bins = rep(rep(factor(levels(cut(0, breaks = seq(0,750,150)))),
                                   length(factor(levels(cut(0, breaks = seq(0,750,150)))))),5),
                       cycle = c(rep(0,25),rep(1,25),rep(2,25),rep(3,25),rep(4,25)))

# create bins in data
treelists_bin <- alldata_spgrp %>%
  mutate(X_bins = cut(X, breaks=seq(0,750,150)),
         Y_bins = cut(Y, breaks=seq(0,750,150)))


#get dominance by stem count (+ functional group and average height)
dom_count <-  treelists_bin %>%
  group_by(cycle,X_bins,Y_bins) %>%
  summarize(avg_height = mean(HT_m),
            func_grp = fmode(func_grp, na.rm=T),
            dom_count = fmode(sci_name,na.rm=T))

#get dominance by biomass
dom_biomass <- treelists_bin %>%
  group_by(cycle,X_bins,Y_bins,sci_name) %>%
  summarize(biomass_sum = sum(AGB_g, na.rm=T)) %>%
  ungroup() %>%
  group_by(cycle,X_bins,Y_bins) %>%
  slice_max(biomass_sum, n=1, with_ties=F) %>%
  rename(dom_biomass=sci_name) %>%
  select(cycle,X_bins,Y_bins,dom_biomass)

#get dominace by BA
dom_BA <- treelists_bin %>%
  group_by(cycle,X_bins,Y_bins,sci_name) %>%
  summarize(DBH_sum = sum(DBH, na.rm=T)) %>%
  ungroup() %>%
  group_by(cycle,X_bins,Y_bins) %>%
  slice_max(DBH_sum, n=1, with_ties=F) %>%
  rename(dom_BA=sci_name) %>%
  select(cycle,X_bins,Y_bins,dom_BA)

#combine them
overlap_dom <- left_join(dom_count,dom_biomass, by = join_by("cycle","X_bins","Y_bins"))
overlap_dom <- left_join(overlap_dom,dom_BA, by=join_by("cycle","X_bins","Y_bins"))


dom_grid <- left_join(all_bins, overlap_dom, by=join_by(cycle,X_bins,Y_bins))

rst_list <- list()
for(i in 0:4){
  BA_mat <- matrix(dom_grid[dom_grid$cycle==i,]$dom_BA,
                   nrow = 5,
                   ncol = 5,
                   byrow = F)
  BA_rst <- flip(rast(xmin = 0, xmax = 750,
                      ymin = 0, ymax = 750,
                      nrow = 5,
                      ncol = 5,
                      vals=BA_mat))
  rst_list[[i+1]] <- BA_rst
}

BA_rasters <- c(rst_list[[1]],rst_list[[2]],rst_list[[3]],rst_list[[4]],rst_list[[5]])

names(BA_rasters) <- c("10 Years\n(Spinup)",
                       "20 Years\n(Fire cycle 1)",
                       "30 Years\n(Fire cycle 2)",
                       "40 Years\n(Fire cycle 3)",
                       "50 Years\n(Fire cycle 4)")
BA_rasters <- BA_rasters %>%
  mutate(`10 Years\n(Spinup)` = factor(`10 Years\n(Spinup)`, 
                             levels = c("Quercus chrysolepis",
                                        "Pseudotsuga menziesii",
                                        "Calocedrus decurrens",
                                        "Pinus ponderosa",
                                        "Pinus contorta",
                                        "Taxus brevifolia")),
         `20 Years\n(Fire cycle 1)` = factor(`20 Years\n(Fire cycle 1)`, 
                             levels = c("Quercus chrysolepis",
                                        "Pseudotsuga menziesii",
                                        "Calocedrus decurrens",
                                        "Pinus ponderosa",
                                        "Pinus contorta",
                                        "Taxus brevifolia")),
         `30 Years\n(Fire cycle 2)` = factor(`30 Years\n(Fire cycle 2)`, 
                             levels = c("Quercus chrysolepis",
                                        "Pseudotsuga menziesii",
                                        "Calocedrus decurrens",
                                        "Pinus ponderosa",
                                        "Pinus contorta",
                                        "Taxus brevifolia")),
         `40 Years\n(Fire cycle 3)` = factor(`40 Years\n(Fire cycle 3)`, 
                             levels = c("Quercus chrysolepis",
                                        "Pseudotsuga menziesii",
                                        "Calocedrus decurrens",
                                        "Pinus ponderosa",
                                        "Pinus contorta",
                                        "Taxus brevifolia")),
         `50 Years\n(Fire cycle 4)` = factor(`50 Years\n(Fire cycle 4)`, 
                             levels = c("Quercus chrysolepis",
                                        "Pseudotsuga menziesii",
                                        "Calocedrus decurrens",
                                        "Pinus ponderosa",
                                        "Pinus contorta",
                                        "Taxus brevifolia")))
ext(BA_rasters) = ext(c(0,750,0,750))


perim1 <- read_perimeters(1)
perim2 <- read_perimeters(2)
perim3 <- read_perimeters(3)
perim4 <- read_perimeters(4)
perimeters <- rbind(perim1,perim2,perim3,perim4)
perimeters <- perimeters %>%
  mutate(lyr = case_when(cycle==1 ~ "10 Years\n(Spinup)",
                         cycle==2 ~ "20 Years\n(Fire cycle 1)",
                           cycle==3 ~ "30 Years\n(Fire cycle 2)",
                           TRUE ~ "40 Years\n(Fire cycle 3)"))
perim_sf <- st_as_sf(perimeters)

## Tree height
rst_list <- list()
for(i in 0:4){
  HT_mat <- matrix(dom_grid[dom_grid$cycle==i,]$avg_height,
                   nrow = 5,
                   ncol = 5,
                   byrow = F)
  HT_rst <- flip(rast(xmin = 0, xmax = 750,
                      ymin = 0, ymax = 750,
                      nrow = 5,
                      ncol = 5,
                      vals=HT_mat))
  rst_list[[i+1]] <- HT_rst
}

height <- c(rst_list[[1]],rst_list[[2]],rst_list[[3]],rst_list[[4]],rst_list[[5]])

names(height) <- c("10 Years\n(Spinup)",
                   "20 Years\n(Fire cycle 1)",
                   "30 Years\n(Fire cycle 2)",
                   "40 Years\n(Fire cycle 3)",
                   "50 Years\n(Fire cycle 4)")
ext(height) = ext(c(0,750,0,750))

dom_species <- ggplot() +
  geom_spatraster(data=BA_rasters) +
  geom_sf_pattern(data=perim_sf, 
                  color = "red", 
                  fill=NA, 
                  pattern = "stripe",
                  pattern_fill = "red",
                  pattern_color = NA,
                  pattern_alpha = 0.7,
                  pattern_size = 0.5) +
  scale_fill_colorblind() +
  scale_y_continuous(breaks = c(0,150,300,450,600,750)) +
  facet_wrap(~lyr, nrow=1) +
  labs(y = "Y (m)",
       fill = "Dominant Species\nby Basal Area") +
  theme(legend.position="bottom",
        legend.text = element_text(face="italic"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank()) +
  guides(fill = guide_legend(nrow=2))

tree_height <- ggplot() +
  geom_spatraster(data=height) +
  geom_sf_pattern(data=perim_sf, 
                  aes(color = factor(IN)), 
                  fill=NA, 
                  pattern = "stripe",
                  pattern_fill = "red",
                  pattern_color = NA,
                  pattern_alpha = 0.7,
                  pattern_size = 0.5) +
  scale_fill_distiller(palette = 2, direction=1, na.value="white") +
  scale_color_manual(values = "red") +
  scale_y_continuous(breaks = c(0,150,300,450,600,750)) +
  scale_x_continuous(breaks = c(0,150,300,450,600,750)) +
  facet_wrap(~lyr, nrow=1) +
  labs(x = "X (m)",
       y = "Y (m)",
       fill = "Average Tree Height (m)",
       color = "Fire Perimeter") +
  theme(legend.position="bottom",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(color = guide_legend(theme = theme(legend.text = element_blank())))

dom_species / tree_height

ggsave("figure_6.png",
       path = here("1.LANDIS-MODEL","figures"), 
       width = col2, 
       height = 5,
       units = "in")

# plot individual trees?
all_postlandis <- alldata_spgrp %>% filter(stage == "postlandis") %>%
  arrange(DBH)

dom_onecell <- all_postlandis %>%
  filter(Y>=300,X>=300) %>%
  filter(Y<450,X<450)

pl_vect <- vect(dom_onecell, geom = c("X","Y"), crs = "") %>%
  rename(Species = sci_name) %>%
  mutate(cycle = factor(cycle, labels = c("10 Years\n(Spinup)",
                                          "20 Years\n(Fire cycle 1)",
                                          "30 Years\n(Fire cycle 2)",
                                          "40 Years\n(Fire cycle 3)",
                                          "50 Years\n(Fire cycle 4)")))
onecell_colors <- colorblind_pal()(8)[c(7,8,3,4,2)]

stems <- ggplot() +
  geom_spatvector(data = pl_vect,
                aes(color = Species, 
                    fill = Species,
                    size = DBH),
                shape = 21,
                alpha = 0.5) +
  scale_size(range = c(0.1,2), 
             limits = c(0,40), 
             breaks = c(10,20,30,40)) +
  scale_color_manual(values = onecell_colors) +
  scale_fill_manual(values = onecell_colors) +
  scale_x_continuous(breaks = c(325,375,425)) +
  scale_y_continuous(breaks = c(325,375,425)) +
  facet_wrap(~cycle, nrow=1) +
  labs(x="X (m)",
       y = "Y (m)",
       size = "DBH (cm)") +
  theme(legend.position = "bottom",
        legend.text = element_text(face = "italic"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=3), nrow=2),
         size = guide_legend(nrow=2))
stems

ggsave("figure_8.png",path = here("1.LANDIS-MODEL","figures"), width = col2, height = 2.8)

dom_othercell <- all_postlandis %>%
  filter(Y>=0,X>=450) %>%
  filter(Y<150,X<600)

other_vect <- vect(dom_othercell, geom = c("X","Y"), crs = "") %>%
  rename(Species = sci_name) %>%
  mutate(cycle = factor(cycle, labels = c("10 Years\n(Spinup)",
                                          "20 Years\n(Fire cycle 1)",
                                          "30 Years\n(Fire cycle 2)",
                                          "40 Years\n(Fire cycle 3)",
                                          "50 Years\n(Fire cycle 4)")))

# take out bigleaf maple since there is barely any and it messes with the aesthetics
other_vect <- other_vect %>% filter(Species != "Acer macrophyllum")
othercell_colors <- colorblind_pal()(8)[c(8,3,4,2,1)]

ggplot() +
  geom_spatvector(data = other_vect,
                  aes(color = Species, 
                      fill = Species,
                      size = DBH),
                  shape = 21,
                  alpha = 0.5) +
  scale_size(range = c(0.1,2)) +
  scale_color_manual(values = othercell_colors) +
  scale_fill_manual(values = othercell_colors) +
  scale_x_continuous(breaks = c(475,525,575)) +
  scale_y_continuous(breaks = c(25,75,125)) +
  facet_wrap(~cycle, nrow=1) +
  labs(x="X (m)",
       y = "Y (m)",
       size = "DBH (cm)") +
  theme(legend.position = "bottom",
        legend.text = element_text(face = "italic"),
        strip.text = element_blank(),
        panel.grid = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=3), nrow=2),
         size = guide_legend(nrow=2))

ggsave("figure_8_alt.png",path = here("1.LANDIS-MODEL","figures"), width = col2, height = 2.8)

## Combining these ^^ in powerpoint

## figure 7
alldata_stems <- alldata_spgrp %>%
  filter(stage == "postlandis") %>%
  mutate(cycle_name = paste("Cycle", cycle)) %>%
  filter(COMMON_NAME %in% c("lodgepole pine",
                            "ponderosa pine",
                            "incense-cedar",
                            "Douglas-fir",
                            "canyon live oak")) %>%
  mutate(sci_name = factor(sci_name, levels = c("Quercus chrysolepis",
                                                "Pseudotsuga menziesii",
                                                "Calocedrus decurrens",
                                                "Pinus ponderosa",
                                                "Pinus contorta"))) %>%
  ggplot() +
  geom_histogram(aes(DBH, fill=sci_name), bins = 12) +
  facet_grid(sci_name ~ cycle_name, scales = "free_y") +
  scale_fill_colorblind() +
  theme_bw() +
  theme(legend.position="none")
alldata_stems



ha = (res(OC)[1]*dim(OC)[1])*(res(OC)[2]*dim(OC)[2])/10000
BA_plot <- alldata_spgrp %>%
  filter(stage == "postlandis") %>%
  filter(COMMON_NAME %in% c("lodgepole pine",
                            "ponderosa pine",
                            "Douglas-fir",
                            "incense-cedar",
                            "canyon live oak")) %>%
  mutate(sci_name = factor(sci_name, levels = c("Quercus chrysolepis",
                                                "Pseudotsuga menziesii",
                                                "Calocedrus decurrens",
                                                "Pinus ponderosa",
                                                "Pinus contorta"))) %>%
  mutate(DBH_bins = cut(DBH, 
                        breaks=c(0,5,10,15,20,25,30,35,40,45,50),
                        labels = c(5,10,15,20,25,30,35,40,45,50))) %>%
  mutate(DBH_bins_num = case_when(DBH_bins==5 ~ 0,
                                  DBH_bins==10 ~ 5,
                                  DBH_bins==15 ~ 10,
                                  DBH_bins==20 ~ 15,
                                  DBH_bins==25 ~ 20,
                                  DBH_bins==30 ~ 25,
                                  DBH_bins==35 ~ 30,
                                  DBH_bins==45 ~ 35,
                                  DBH_bins==45 ~ 40,
                                  DBH_bins==50 ~ 45,
                                  TRUE ~ 50)) %>%
  mutate(BA = (pi*(DBH/2)^2)/10000) %>%
  mutate(cycle = factor(cycle, labels = c("10 Years\n(Spinup)",
                                          "20 Years\n(Fire cycle 1)",
                                          "30 Years\n(Fire cycle 2)",
                                          "40 Years\n(Fire cycle 3)",
                                          "50 Years\n(Fire cycle 4)"))) %>%
  group_by(cycle, sci_name, DBH_bins) %>%
  summarize(BA_sum = sum(BA)/ha,
            DBH_bins = mean(DBH_bins_num)) %>%
  ggplot() +
  geom_bar(stat="identity", 
           aes(DBH_bins, BA_sum, fill=sci_name),
           just=0,
           width = 5) +
  facet_grid(sci_name ~ cycle) +
  scale_fill_colorblind() +
  labs(x="Diameter at Breast Height (cm)",
       y = "Basal Area (sq. m per ha)") +
  theme(legend.position="none")
BA_plot

ggsave("figure_7.png",path = here("1.LANDIS-MODEL","figures"), width = col2, height = 7.7)

#############
## how many trees died in each fire?
all_colors = c(colorblind_pal()(8),"tan")
trees_killed <- alldata_spgrp %>%
  mutate(fire_ID = if_else(stage=="postlandis",cycle+1,cycle)) %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus chrysolepis",
                                            "Pseudotsuga menziesii",
                                            "Calocedrus decurrens",
                                            "Pinus ponderosa",
                                            "Pinus contorta",
                                            "Taxus brevifolia",
                                            "Abies grandis",
                                            "Arbutus menziesii",
                                            "Other hardwood"))) %>%
  group_by(stage, fire_ID, sci_name_other) %>%
  summarize(n_trees = n()) %>%
  pivot_wider(names_from = stage,
              values_from = n_trees) %>%
  mutate(n_dead = (postlandis - postfire)) %>%
  filter(fire_ID != 5) %>%
  ggplot() +
  geom_bar(aes(fire_ID, n_dead, fill = sci_name_other),
              stat = "identity",
           position = "stack") +
  scale_fill_manual(values = all_colors) +
  labs(x = "Fire cycle",
       y = "Trees killed",
       fill = "Species") +
  theme(legend.text = element_text(face = "italic"),
        axis.text.y = element_text(angle = 45))
trees_killed

ggsave("trees_killed.png",path = here("1.LANDIS-MODEL","figures"), width = col1, height = 3)

# biomass loss
biomass_lost <- alldata_spgrp %>%
  mutate(fire_ID = if_else(stage=="postlandis",cycle+1,cycle)) %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus chrysolepis",
                                            "Pseudotsuga menziesii",
                                            "Calocedrus decurrens",
                                            "Pinus ponderosa",
                                            "Pinus contorta",
                                            "Taxus brevifolia",
                                            "Abies grandis",
                                            "Arbutus menziesii",
                                            "Other hardwood"))) %>%
  group_by(stage, fire_ID, sci_name_other) %>%
  summarize(biomass_Mg = sum(AGB_g)/1000000) %>%
  pivot_wider(names_from = stage,
              values_from = biomass_Mg) %>%
  mutate(biomass_loss = postlandis - postfire) %>%
  filter(fire_ID != 5) %>%
  ggplot() +
  geom_bar(aes(fire_ID, biomass_loss, fill = sci_name_other),
           stat = "identity",
           position = "stack") +
  scale_fill_manual(values = all_colors) +
  labs(x = "Fire cycle",
       y = "Aboveground biomass loss (Mg)",
       fill = "Species") +
  theme(legend.text = element_text(face = "italic"),
        axis.text.y = element_text(angle = 45))
biomass_lost

ggsave("biomass_loss.png",path = here("1.LANDIS-MODEL","figures"), width = col1, height = 3)

# proportion total biomass
total_biomass <- alldata_spgrp %>%
  mutate(fire_ID = if_else(stage=="postlandis",cycle+1,cycle)) %>%
  group_by(stage, fire_ID) %>%
  summarize(total_biomass_Mg = sum(AGB_g)/1000000)

prop_biomass <- alldata_spgrp %>%
  mutate(fire_ID = if_else(stage=="postlandis",cycle+1,cycle)) %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus chrysolepis",
                                            "Pseudotsuga menziesii",
                                            "Calocedrus decurrens",
                                            "Pinus ponderosa",
                                            "Pinus contorta",
                                            "Taxus brevifolia",
                                            "Abies grandis",
                                            "Arbutus menziesii",
                                            "Other hardwood"))) %>%
  group_by(stage, fire_ID, sci_name_other) %>%
  summarize(biomass_Mg = sum(AGB_g)/1000000) %>%
  left_join(total_biomass, by = join_by("fire_ID","stage")) %>%
  mutate(prop_biomass = (biomass_Mg/total_biomass_Mg)*100,
         fire_ID = fire_ID*10) %>%
  filter(stage=="postlandis") %>%
  ggplot() +
  geom_bar(aes(fire_ID, prop_biomass, fill = sci_name_other),
           stat = "identity",
           position = "stack") +
  scale_fill_manual(values = all_colors) +
  labs(x = "Years",
       y = "Total biomass (%)",
       fill = "Species") +
  theme(legend.text = element_text(face = "italic"))
prop_biomass

ggsave("proportion_biomass.png",path = here("1.LANDIS-MODEL","figures"), width = col1, height = 3)


# sizes of trees that died

pl_spgrp1 <- alldata_spgrp %>%
  filter(stage == "postlandis") %>%
  mutate(fire_ID = cycle+1)
pl_spgrp2 <- pl_spgrp1 %>%
  select(X,Y,fire_ID) %>%
  mutate(prefire = 1)
pf_spgrp <- alldata_spgrp %>% 
  filter(stage == "postfire") %>%
  rename(fire_ID = cycle) %>%
  select(X,Y,fire_ID) %>%
  mutate(postfire = 1)
  
master = left_join(pl_spgrp2, pf_spgrp)
dead_trees <- master %>% 
  filter(is.na(postfire)) %>% 
  select(-prefire,-postfire) %>%
  left_join(pl_spgrp1) %>%
  filter(fire_ID != 5)

dead_trees %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = str_replace(sci_name_other, " ", "\n")) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus\nchrysolepis",
                                            "Pseudotsuga\nmenziesii",
                                            "Calocedrus\ndecurrens",
                                            "Pinus\nponderosa",
                                            "Pinus\ncontorta",
                                            "Taxus\nbrevifolia",
                                            "Abies\ngrandis",
                                            "Arbutus\nmenziesii",
                                            "Other\nhardwood"))) %>%
  mutate(fire_ID = factor(fire_ID, labels = c("Fire cycle 1",
                                              "Fire cycle 2",
                                              "Fire cycle 3",
                                              "Fire cycle 4"))) %>%
  ggplot() +
  geom_histogram(aes(DBH,
                     after_stat(count)/tapply(after_stat(count), after_stat(PANEL), sum)[after_stat(PANEL)], 
                     fill = sci_name_other),
                 binwidth = 5) +
  scale_fill_manual(values = all_colors) +
  facet_grid(sci_name_other ~ fire_ID, scales = "free") +
  scale_y_continuous(limits = c(0,1), 
                     labels = percent_format()) +
  labs(x = "DBH (cm)",
       y = "Trees killed in each fire",
       fill = "Species") +
  theme(legend.position = "none",
        strip.text.y = element_text(face = "italic"))

ggsave("dbh_killed.png", path = here("1.LANDIS-MODEL","figures"),width=col2, height = 9)

# stacked biomass through time
alldata_spgrp %>%
  mutate(fire_ID = if_else(stage=="postlandis",cycle+1,cycle)) %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus chrysolepis",
                                            "Pseudotsuga menziesii",
                                            "Calocedrus decurrens",
                                            "Pinus ponderosa",
                                            "Pinus contorta",
                                            "Taxus brevifolia",
                                            "Abies grandis",
                                            "Arbutus menziesii",
                                            "Other hardwood"))) %>%
  group_by(stage, fire_ID, sci_name_other) %>%
  summarize(biomass_Mg = sum(AGB_g)/1000000) %>%
  mutate(timestep = case_when(fire_ID==1 & stage=="postlandis" ~ "10 Years\nPrefire",
                              fire_ID==1 & stage=="postfire" ~ "10 Years\nPostfire",
                              fire_ID==2 & stage=="postlandis" ~ "20 Years\nPrefire",
                              fire_ID==2 & stage=="postfire" ~ "20 Years\nPostfire",
                              fire_ID==3 & stage=="postlandis" ~ "30 Years\nPrefire",
                              fire_ID==3 & stage=="postfire" ~ "30 Years\nPostfire",
                              fire_ID==4 & stage=="postlandis" ~ "40 Years\nPrefire",
                              fire_ID==4 & stage=="postfire" ~ "40 Years\nPostfire",
                              TRUE ~ "50 Years")) %>%
  mutate(timestep = factor(timestep,
                           levels = c("10 Years\nPrefire",
                                      "10 Years\nPostfire",
                                      "20 Years\nPrefire",
                                      "20 Years\nPostfire",
                                      "30 Years\nPrefire",
                                      "30 Years\nPostfire",
                                      "40 Years\nPrefire",
                                      "40 Years\nPostfire",
                                      "50 Years"))) %>%
  ggplot() +
  geom_bar(aes(timestep, biomass_Mg, fill = sci_name_other),
           stat = "identity",
           position = "stack") +
  facet_wrap(~ timestep, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = all_colors) +
  scale_y_continuous(limits = c(0,8000),expand=c(0,0)) +
  labs(x = "Years",
       y = "Biomass (Mg)",
       fill = "Species") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        panel.spacing = grid::unit(-0.5, "lines"),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_blank())

ggsave("biomass_thru_time.png",path = here("1.LANDIS-MODEL","figures"), width = 8, height = 3)

# biomass change after fires
biomass_change <- alldata_spgrp %>%
  filter(stage=="postfire") %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus chrysolepis",
                                            "Pseudotsuga menziesii",
                                            "Calocedrus decurrens",
                                            "Pinus ponderosa",
                                            "Pinus contorta",
                                            "Taxus brevifolia",
                                            "Abies grandis",
                                            "Arbutus menziesii",
                                            "Other hardwood"))) %>%
  group_by(cycle, sci_name_other) %>%
  summarize(biomass_Mg = sum(AGB_g)/1000000) %>%
  ggplot() +
  geom_bar(aes(cycle, biomass_Mg, fill = sci_name_other),
           stat = "identity",
           position = "stack") +
  scale_fill_manual(values = all_colors) +
  labs(x = "Fire cycle",
       y = "Aboveground biomass (Mg)",
       fill = "Species") +
  theme(legend.text = element_text(face = "italic"),
        axis.text.y = element_text(angle = 45))

trees_killed + biomass_lost + biomass_change +
  plot_layout(nrow = 1,
              guides = "collect") &
  theme(legend.background = element_rect(fill = NA))
ggsave("fire_effects.png",path = here("1.LANDIS-MODEL","figures"), width = col2, height = 3)
  
## quicfire
qf <- read.csv(here("1.LANDIS-MODEL","QF_outputs","qf_consumption.csv"))

consumption <- qf %>%
  pivot_longer(cols = 2:9,
               names_to = "var",
               values_to = "val") %>%
  mutate(strata = case_when(var %in% c("surf_cons","surf_cons_pct") ~ "Surface\n(<1 m)",
                            var %in% c("canopy_cons","can_cons_pct") ~ "Canopy\n(>1 m)",
                            var %in% c("midstory_cons","mid_cons_pct") ~ "Ladder fuels\n(1-7 m)",
                            var %in% c("overstory_cons","over_cons_pct") ~ "Overstory\n(>7 m)"),
         unit = case_when(var %in% c("surf_cons",
                                     "canopy_cons",
                                     "midstory_cons",
                                     "overstory_cons") ~ "Total (Mg)",
                          TRUE ~ "Percent")) %>%
  mutate(val = if_else(unit=="Percent",val,val/1000),
         unit = factor(unit, levels = c("Total (Mg)","Percent")),
         strata = factor(strata, levels = c("Surface\n(<1 m)",
                                            "Canopy\n(>1 m)",
                                            "Ladder fuels\n(1-7 m)",
                                            "Overstory\n(>7 m)"))) %>%
  filter(strata != "Canopy\n(>1 m)")

tot_cons <- consumption %>%
  filter(unit == "Total (Mg)") %>%
  ggplot() +
  geom_line(aes(cycle,val)) +
  geom_point(aes(cycle,val), alpha = 0.75) +
  facet_wrap(~strata) +
  scale_y_continuous(limits = c(0,72),
                     expand = c(0,0)) +
  labs(x = "Fire cycle",
       y = "Total (Mg)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor.x = element_blank())

pct_cons <- consumption %>%
  filter(unit == "Percent") %>%
  ggplot() +
  geom_line(aes(cycle,val)) +
  geom_point(aes(cycle,val), alpha = 0.75) +
  facet_wrap(~strata) +
  scale_y_continuous(limits = c(0,1.02),
                     expand = c(0,0),
                     labels = percent_format()) +
  labs(x = "Fire cycle",
       y = "Percent") +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.minor.x = element_blank())

cons <- tot_cons / pct_cons + plot_layout(axis = "collect")
full_cons = patchwork::patchworkGrob(cons)
full_cons <- gridExtra::grid.arrange(full_cons, left = "Fine fuel consumption")

ggsave("fuel_consumption.png",
       plot = full_cons, 
       path = here("1.LANDIS-MODEL","figures"), 
       width = col1, 
       height = 3.2)

###################

topleft <- alldata_spgrp %>%
  filter(X < 300, Y > 450) %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus chrysolepis",
                                            "Pseudotsuga menziesii",
                                            "Calocedrus decurrens",
                                            "Pinus ponderosa",
                                            "Pinus contorta",
                                            "Taxus brevifolia",
                                            "Abies grandis",
                                            "Arbutus menziesii",
                                            "Other hardwood")))
topleft %>%
  filter(stage == "postlandis") %>%
  arrange(DBH) %>%
  ggplot() +
  geom_point(aes(X,Y,color = sci_name_other, size=DBH), alpha = 0.5) +
  scale_color_manual(values=colorblind_pal()(8)[c(2,3,4,5)]) +
  scale_size(range=c(0.1,2)) +
  facet_wrap(~cycle, nrow=1) +
  coord_equal() +
  theme(legend.position = "bottom")

alldata_spgrp %>%
  filter(stage == "postlandis") %>%
  filter(cycle%in%c(0,4)) %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus chrysolepis",
                                            "Pseudotsuga menziesii",
                                            "Calocedrus decurrens",
                                            "Pinus ponderosa",
                                            "Pinus contorta",
                                            "Taxus brevifolia",
                                            "Abies grandis",
                                            "Arbutus menziesii",
                                            "Other hardwood"))) %>%
  arrange(DBH) %>%
  mutate(cycle = factor(cycle, labels = c("Initial conditions", "After 50 yr. DRM"))) %>%
  ggplot() +
  geom_point(aes(X,Y,color = sci_name_other, size=DBH), alpha = 0.5) +
  scale_color_manual(values=all_colors) +
  scale_size(range=c(0.005,1.2)) +
  facet_wrap(~cycle, nrow=1) +
  labs(x = "X (m)",
       y = "Y (m)",
       size = "DBH (cm)",
       color = "Species") +
  coord_equal() +
  theme(legend.text = element_text(face="italic"),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3),
         size = guide_legend(nrow = 3))

ggsave("before_after.png",
       path = here("1.LANDIS-MODEL","figures"), 
       width = col2, 
       height = 5)


##########
# for each size class, show total number of trees, split up by which survived and which died (separate panels)
# split this up by species

pl_spgrp1 <- alldata_spgrp %>%
  filter(stage == "postlandis") %>%
  mutate(fire_ID = cycle+1)
pl_spgrp2 <- pl_spgrp1 %>%
  select(X,Y,fire_ID) %>%
  mutate(prefire = 1)
pf_spgrp <- alldata_spgrp %>% 
  filter(stage == "postfire") %>%
  rename(fire_ID = cycle) %>%
  select(X,Y,fire_ID) %>%
  mutate(postfire = 1)

master = left_join(pl_spgrp2, pf_spgrp)

dead_live <- master %>%
  mutate(live_dead = if_else(is.na(postfire),"Killed","Survived")) %>%
  select(-prefire,-postfire) %>%
  left_join(pl_spgrp1) %>%
  filter(fire_ID != 5)

strip_colors = all_colors
strip_colors[1] <- "gray80"
strip <- strip_themed(background_y = elem_list_rect(fill = strip_colors))

dead_live %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = str_replace(sci_name_other, " ", "\n")) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus\nchrysolepis",
                                            "Pseudotsuga\nmenziesii",
                                            "Calocedrus\ndecurrens",
                                            "Pinus\nponderosa",
                                            "Pinus\ncontorta",
                                            "Taxus\nbrevifolia",
                                            "Abies\ngrandis",
                                            "Arbutus\nmenziesii",
                                            "Other\nhardwood"))) %>%
  mutate(fire_ID = factor(fire_ID, 
                          levels = c(1,2,3,4),
                          labels = c("Fire 1",
                                     "Fire 2",
                                     "Fire 3",
                                     "Fire 4"))) %>%
  mutate(DBH_bins = cut(DBH, 
                        breaks = seq(0, ceiling(max(DBH)), 
                                     by = 5), 
                        right = FALSE)) %>%
  filter(!is.na(DBH_bins)) %>%
  mutate(bin_start = as.numeric(sub("\\[([-+]?[0-9]*\\.?[0-9]+),.*", "\\1", DBH_bins))) %>%
  group_by(fire_ID,sci_name_other,bin_start) %>%
  summarize(
    Killed = sum(live_dead == "Killed"),
    Survived = sum(live_dead == "Survived"),
    Count = n()
  ) %>%
  pivot_longer(cols = c("Killed","Survived"),
               names_to = "killed_survived",
               values_to = "ks_pct") %>%
  mutate(killed_survived = factor(killed_survived, levels = c("Survived","Killed"))) %>%
  ggplot() +
  geom_col(aes(bin_start, ks_pct, fill = killed_survived), 
           width = 5,
           position = "stack",
           just = 0) +
  facet_grid2(sci_name_other~fire_ID, strip = strip, scales = "free_y") +
  scale_fill_manual(values = c("gray80","gray30")) +
  labs(x = "DBH (cm)",
       y = "Proportion of stems killed or survived") +
  theme(legend.title = element_blank())

ggsave("killed_size_class.png",path = here("1.LANDIS-MODEL","figures"), width = col2, height = 7.7)

dead_live %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = str_replace(sci_name_other, " ", "\n")) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus\nchrysolepis",
                                            "Pseudotsuga\nmenziesii",
                                            "Calocedrus\ndecurrens",
                                            "Pinus\nponderosa",
                                            "Pinus\ncontorta",
                                            "Taxus\nbrevifolia",
                                            "Abies\ngrandis",
                                            "Arbutus\nmenziesii",
                                            "Other\nhardwood"))) %>%
  mutate(fire_ID = factor(fire_ID, 
                          levels = c(1,2,3,4),
                          labels = c("Fire 1",
                                     "Fire 2",
                                     "Fire 3",
                                     "Fire 4"))) %>%
  mutate(DBH_bins = cut(DBH, 
                        breaks = seq(0, ceiling(max(DBH)), 
                                     by = 5), 
                        right = FALSE)) %>%
  filter(!is.na(DBH_bins)) %>%
  mutate(bin_start = as.numeric(sub("\\[([-+]?[0-9]*\\.?[0-9]+),.*", "\\1", DBH_bins))) %>%
  group_by(fire_ID,bin_start) %>%
  summarize(
    Killed = mean(live_dead == "Killed"),
    Survived = mean(live_dead == "Survived"),
    Count = n()
  ) %>%
  ggplot() +
  geom_col(aes(bin_start, Survived),
           fill = "gray30",
           width = 5,
           position = "stack",
           just = 0) +
  facet_wrap(~fire_ID, nrow=1) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "DBH class (cm)",
       y = "Proportion of stems killed or survived") +
  theme(legend.title = element_blank())


###########
cons <- read.csv(here("1.LANDIS-MODEL","QF_outputs","consumption_rate.csv"))
cons$timestep <- row.names(cons)
cons <- cons %>%
  mutate(time = as.numeric(timestep)*30) %>%
  pivot_longer(cols = 1:12,
               names_to = "strata_cycle",
               values_to = "cons_rate")

cons$strata = NA
cons$cycle = NA
for(i in 1:nrow(cons)){
  cons$strata[i] <- str_split(cons$strata_cycle[i],"_")[[1]][1]
  cons$cycle[i] <- as.numeric(str_split(cons$strata_cycle[i],"_")[[1]][2])
}

cons_colors <- solarized_pal()(7)[c(7,5,6,3)]

cons %>%
  select(-strata_cycle) %>%
  filter(time < 1500) %>%
  mutate(strata = factor(strata, 
                        levels = c("surf","mid","over"),
                        labels = c("Surface (<1 m)","Midstory (1-7 m)","Overstory (>7 m)")),
         cycle = factor(cycle,
                        levels = c(1,2,3,4),
                         labels = c("Fire 1","Fire 2","Fire 3","Fire 4"))) %>%
  ggplot() +
  geom_line(aes(time,cons_rate,color=cycle)) +
  facet_wrap(~strata, nrow=1) +
  scale_color_manual(values = cons_colors) +
  scale_x_continuous(breaks = c(0,300,600,900,1200)) +
  labs(x="Time (s)",
       y="Consumption Rate (kg/s)",
       linetype = "Fuel strata") +
  theme(legend.title = element_blank())

ggsave("consumption_rate.png",
       path = here("1.LANDIS-MODEL","figures"), 
       width = col2, 
       height = 2.5)

########
# highlight loss of ladder fuels

alldata_spgrp %>%
  filter(stage == "postfire") %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus chrysolepis",
                                            "Pseudotsuga menziesii",
                                            "Calocedrus decurrens",
                                            "Pinus ponderosa",
                                            "Pinus contorta",
                                            "Taxus brevifolia",
                                            "Abies grandis",
                                            "Arbutus menziesii",
                                            "Other hardwood"))) %>%
  mutate(fire_ID = factor(cycle, 
                          levels = c(1,2,3,4),
                          labels = c("Fire 1",
                                     "Fire 2",
                                     "Fire 3",
                                     "Fire 4"))) %>%
  mutate(canopy_class = if_else(HTLC_m < 7, "Ladder fuel (CBH < 7 m)", "Canopy trees (CBH > 7 m)")) %>%
  mutate(canopy_class = factor(canopy_class, levels = c("Ladder fuel (CBH < 7 m)", "Canopy trees (CBH > 7 m)"))) %>%
  group_by(fire_ID, sci_name_other, canopy_class) %>%
  summarize(n_trees = n()) %>%
  ggplot() +
  geom_bar(aes(fire_ID,n_trees,fill=sci_name_other),
           stat = "identity", 
           position = "stack") +
  facet_wrap(~canopy_class) +
  scale_fill_manual(values = all_colors) +
  labs(x = "Fire",
       y = "Trees remaining after fire",
       fill = "Species") +
  theme(legend.text = element_text(face = "italic"),
        axis.title.x = element_blank())

ggsave("ladder_fuels.png",
       path = here("1.LANDIS-MODEL","figures"), 
       width = col2, 
       height = 3)

##
pl_spgrp1 <- alldata_spgrp %>%
  filter(stage == "postlandis") %>%
  mutate(fire_ID = cycle+1)
pl_spgrp2 <- pl_spgrp1 %>%
  select(X,Y,fire_ID) %>%
  mutate(prefire = 1)
pf_spgrp <- alldata_spgrp %>% 
  filter(stage == "postfire") %>%
  rename(fire_ID = cycle) %>%
  select(X,Y,fire_ID) %>%
  mutate(postfire = 1)

master = left_join(pl_spgrp2, pf_spgrp)

dead_live <- master %>%
  mutate(live_dead = if_else(is.na(postfire),"Killed","Survived")) %>%
  select(-prefire,-postfire) %>%
  left_join(pl_spgrp1) %>%
  filter(fire_ID != 5)

strip_colors = all_colors
strip_colors[1] <- "gray80"
strip <- strip_themed(background_y = elem_list_rect(fill = strip_colors))

dead_live %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = str_replace(sci_name_other, " ", "\n")) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus\nchrysolepis",
                                            "Pseudotsuga\nmenziesii",
                                            "Calocedrus\ndecurrens",
                                            "Pinus\nponderosa",
                                            "Pinus\ncontorta",
                                            "Taxus\nbrevifolia",
                                            "Abies\ngrandis",
                                            "Arbutus\nmenziesii",
                                            "Other\nhardwood"))) %>%
  mutate(fire_ID = factor(fire_ID, 
                          levels = c(1,2,3,4),
                          labels = c("Fire 1",
                                     "Fire 2",
                                     "Fire 3",
                                     "Fire 4"))) %>%
  mutate(canopy_class = if_else(HTLC_m < 7, "Ladder fuel (CBH < 7 m)", "Canopy trees (CBH > 7 m)")) %>%
  mutate(canopy_class = factor(canopy_class, levels = c("Ladder fuel (CBH < 7 m)", "Canopy trees (CBH > 7 m)"))) %>%
  group_by(fire_ID,sci_name_other,canopy_class) %>%
  summarize(
    Killed = mean(live_dead == "Killed"),
    Survived = mean(live_dead == "Survived"),
    Count = n()
  ) %>%
  ggplot() +
  geom_bar(aes(canopy_class, Survived, fill = sci_name_other), 
           stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(values = all_colors) +
  scale_y_continuous(limits = c(0,1),
                     labels = percent_format()) +
  labs(x = "Canopy class",
       y = "Proportion of trees that survived") +
  theme(legend.title = element_blank())

########
# another Adam request

pl_spgrp1 <- alldata_spgrp %>%
  filter(stage == "postlandis") %>%
  mutate(fire_ID = cycle+1)
pl_spgrp2 <- pl_spgrp1 %>%
  select(X,Y,fire_ID) %>%
  mutate(prefire = 1)
pf_spgrp <- alldata_spgrp %>% 
  filter(stage == "postfire") %>%
  rename(fire_ID = cycle) %>%
  select(X,Y,fire_ID) %>%
  mutate(postfire = 1)

master = left_join(pl_spgrp2, pf_spgrp)

dead_live <- master %>%
  mutate(live_dead = if_else(is.na(postfire),"Killed","Survived")) %>%
  select(-prefire,-postfire) %>%
  left_join(pl_spgrp1) %>%
  filter(fire_ID != 5)

strip_colors2 = all_colors[2:5]
strip <- strip_themed(background_x = elem_list_rect(fill = strip_colors2))

dead_live %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = str_replace(sci_name_other, " ", "\n")) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus\nchrysolepis",
                                            "Pseudotsuga\nmenziesii",
                                            "Calocedrus\ndecurrens",
                                            "Pinus\nponderosa",
                                            "Pinus\ncontorta",
                                            "Taxus\nbrevifolia",
                                            "Abies\ngrandis",
                                            "Arbutus\nmenziesii",
                                            "Other\nhardwood"))) %>%
  # mutate(DBH_bins = cut(DBH,
  #                       breaks = c(0,10,30,ceiling(max(DBH))),
  #                       right = FALSE)) %>%
  mutate(crown_class = if_else(HTLC_m < 7, "Ladder fuels\n(CBH < 7m)", "Canopy trees\n(CBH > 7 m)")) %>%
  mutate(crown_class = factor(crown_class, levels = c("Ladder fuels\n(CBH < 7m)", "Canopy trees\n(CBH > 7 m)"))) %>%
  # filter(!is.na(DBH_bins)) %>%
  filter(sci_name_other %in% c("Pseudotsuga\nmenziesii",
                               "Calocedrus\ndecurrens",
                               "Pinus\nponderosa",
                               "Pinus\ncontorta")) %>%
  group_by(fire_ID,sci_name_other,crown_class) %>%
  summarize(
    Killed = sum(live_dead == "Killed"),
    Survived = sum(live_dead == "Survived"),
    Count = n()
  ) %>%
  ggplot() +
  geom_line(aes(fire_ID,Killed,color=crown_class)) +
  facet_wrap2(~sci_name_other, strip = strip, nrow = 1) +
  # scale_y_continuous(labels = percent_format()) +
  labs(x = "Fire",
       # y = "Proportion of stems killed",
       y = "Number of stems killed",
       color = "Canopy class")

ggsave("killed_size_class_lines_cc_stems.png",path = here("1.LANDIS-MODEL","figures"), width = col2, height = 3)

## ###
# proportion biomass vs size class over time

alldata_spgrp %>%
  filter(stage=="postfire") %>%
  mutate(sci_name_other = if_else(sci_name %in% c("Acer macrophyllum",
                                                  "Alnus rubra",
                                                  "Notholithocarpus densiflorus",
                                                  "Quercus kelloggii"),
                                  "Other hardwood",
                                  sci_name)) %>%
  mutate(sci_name_other = str_replace(sci_name_other, " ", "\n")) %>%
  mutate(sci_name_other = factor(sci_name_other,
                                 levels = c("Quercus\nchrysolepis",
                                            "Pseudotsuga\nmenziesii",
                                            "Calocedrus\ndecurrens",
                                            "Pinus\nponderosa",
                                            "Pinus\ncontorta",
                                            "Taxus\nbrevifolia",
                                            "Abies\ngrandis",
                                            "Arbutus\nmenziesii",
                                            "Other\nhardwood"))) %>%
  # mutate(DBH_bins = cut(DBH,
  #                       breaks = c(0,10,30,ceiling(max(DBH))),
  #                       right = FALSE)) %>%
  mutate(crown_class = if_else(HTLC_m < 7, "Ladder fuels\n(CBH < 7m)", "Canopy trees\n(CBH > 7 m)")) %>%
  mutate(crown_class = factor(crown_class, levels = c("Ladder fuels\n(CBH < 7m)", "Canopy trees\n(CBH > 7 m)"))) %>%
  # filter(!is.na(DBH_bins)) %>%
  filter(sci_name_other %in% c("Pseudotsuga\nmenziesii",
                               "Calocedrus\ndecurrens",
                               "Pinus\nponderosa",
                               "Pinus\ncontorta")) %>%
  group_by(cycle,crown_class) %>%
  # group_by(cycle,DBH_bins) %>%
  summarize(AGB_sum = sum(AGB_g)/1000000) %>%
  # ungroup() %>%
  # group_by(cycle) %>%
  # mutate(AGB_pct = AGB_sum/sum(AGB_sum)) %>%
  ggplot() +
  geom_bar(
    aes(cycle,AGB_sum,fill=crown_class),
    # aes(cycle,AGB_sum,fill=DBH_bins),
           stat = "identity",
           position = "stack") +
  # facet_wrap2(~sci_name_other, strip = strip, nrow = 1) +
  # scale_y_continuous(labels = percent_format()) +
  labs(x = "Fire",
       y = "Aboveground biomass (Mg)",
       fill = "Canopy class"
       # fill = "DBH class (cm)"
       )

ggsave("total_agb_cc.png",path = here("1.LANDIS-MODEL","figures"), width = col1, height = 3)

############
# avg tree height before fire
# avg tree height after fire
# avg cbh before fire
# avg cbh after fire

ht_cbh <- alldata_spgrp %>%
  mutate(fire_ID = if_else(stage=="postlandis", cycle+1, cycle)) %>%
  filter(fire_ID != 5) %>%
  mutate(fire_ID = factor(fire_ID, 
                          labels = c("Fire 1","Fire 2", "Fire 3", "Fire 4")),
         stage = factor(stage, 
                        levels = c("postlandis","postfire"),
                        labels = c("Prefire","Postfire"))) %>%
  group_by(stage,fire_ID) %>%
  summarize(HT_avg = mean(HT_m),
            HTLC_avg = mean(HTLC_m),
            HT_se = sd(HT_m)/sqrt(n()),
            HTLC_se = sd(HTLC_m)/sqrt(n()),
            HT_sd = sd(HT_m),
            HTLC_sd = sd(HTLC_m),
            TPH = n()/((750*750)/(100*100))) %>%
  pivot_longer(cols = c("HT_avg","HT_se","HTLC_avg","HTLC_se","HT_sd","HTLC_sd"), 
               names_to = c("var", "stat"), 
               names_sep = "_", 
               values_to = "val") %>%
  pivot_wider(names_from = stat, values_from = val)

p1 <- ht_cbh %>%
  filter(var=="HT") %>%
  ggplot() +
  geom_bar(aes(x=fire_ID, 
                 y=avg, 
                 fill = fire_ID), 
           stat = "identity",
           width = 0.8) +
  # geom_errorbar(aes(x=fire_ID, 
  #                   ymin = avg - se, 
  #                   ymax = avg + se), 
  #               width = 0.3) +
  scale_x_discrete(limits = rev(levels(ht_cbh$fire_ID))) +
  scale_fill_manual(values = cons_colors) +
  labs(x = "Fire",
       y = "Average\ntree height (m)") +
  facet_wrap(~stage, ncol=1, strip.position = "right") +
  coord_flip() +
  theme(axis.title.y= element_blank(),
        legend.position = "none",
        strip.text = element_blank())

p2 <- ht_cbh %>%
  filter(var=="HTLC") %>%
  ggplot() +
  geom_bar(aes(x=fire_ID,
                 y=avg, 
                 fill = fire_ID),
           stat = "identity",
           width = 0.8) +
  # geom_errorbar(aes(x=fire_ID, 
  #                   ymin = avg - se, 
  #                   ymax = avg + se), 
  #               width = 0.3) +
  scale_x_discrete(limits = rev(levels(ht_cbh$fire_ID))) +
  scale_fill_manual(values = cons_colors) +
  labs(x = "Fire",
       y = "Average canopy\nbase height (m)") +
  facet_wrap(~stage, ncol=1, strip.position = "right") +
  coord_flip() +
  theme(axis.title.y= element_blank(),
        legend.position = "none",
        strip.text = element_blank())

p3 <- ht_cbh %>%
  group_by(stage, fire_ID) %>%
  summarize(TPH = mean(TPH)) %>%
  ggplot() +
  geom_bar(aes(x=fire_ID,
               y=TPH, 
               fill = fire_ID),
           stat = "identity",
           width = 0.8) +
  scale_x_discrete(limits = rev(levels(ht_cbh$fire_ID))) +
  scale_fill_manual(values = cons_colors) +
  labs(x = "Fire",
       y = "Trees per ha") +
  facet_wrap(~stage, ncol=1, strip.position = "right") +
  coord_flip() +
  theme(axis.title.y= element_blank(),
        legend.position = "none")

p1 + p2 + plot_layout(nrow = 1, axes = "collect", axis_titles = "collect")

ggsave("ht_cbh.png",path = here("1.LANDIS-MODEL","figures"), width = col1, height = 3)

p1 + p2 + p3 + plot_layout(nrow = 1, axes = "collect", axis_titles = "collect")

ggsave("ht_cbh_tph.png",path = here("1.LANDIS-MODEL","figures"), width = col2, height = 3)

#########
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
  theme(legend.title = element_blank())

ggsave("windspeeds.png",path = here("1.LANDIS-MODEL","figures"), width = col1, height = 3)
