library(here)
library(tidyverse)
library(terra)
library(ggthemes)
library(tidyterra)
library(collapse)
library(patchwork)
library(sf)
library(gridExtra)

cycle_rst <- function(cyc, stage, rast_df, factor_df, alldata_df){
  t_rst_df <- left_join(rast_df, factor_df[factor_df$cycle==cyc & factor_df$stage %in% stage,], by=join_by("MapCode"))
  count_labels <- unique(alldata_df[alldata_df$cycle==cyc 
                                    & factor_df$stage %in% stage
                                    & alldata_df$MapCode %in% t_rst_df$MapCode,]$dom_count)
  biomass_labels <- unique(alldata_df[alldata_df$cycle==cyc 
                                      & factor_df$stage %in% stage
                                      & alldata_df$MapCode %in% t_rst_df$MapCode,]$dom_biomass)
  BA_labels <- unique(alldata_df[alldata_df$cycle==cyc 
                                 & factor_df$stage %in% stage 
                                 & alldata_df$MapCode %in% t_rst_df$MapCode,]$dom_BA)
  t_rst <- rast(t_rst_df) %>%
    mutate(
      func_grp = factor(
        func_grp,
        labels = c("Conifer","Hardwood")
      ),
      dom_count = factor(
        dom_count, 
        labels = count_labels
      ),
      dom_biomass = factor(
        dom_biomass, 
        labels = biomass_labels
      ),
      dom_BA = factor(
        dom_BA, 
        labels = BA_labels
      )
    )
  return(t_rst)
}

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
  mutate(X_bins = cut(X, breaks=seq(0,1500,150)),
         Y_bins = cut(Y, breaks=seq(0,1800,150)))


#get dominance by stem count (+ functional group and average height)
dom_count <-  treelists_bin %>%
  group_by(cycle,X_bins,Y_bins) %>%
  summarize(avg_height = mean(HT_m),
            func_grp = fmode(func_grp, na.rm=T),
            dom_count = fmode(sci_name,na.rm=T))

#get dominace by biomass
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
ext(basal_area) = ext(c(0,750,0,750))


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
  theme_bw() +
  theme(legend.position="bottom",
        legend.text = element_text(face="italic"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank()) +
  guides(fill = guide_legend(nrow=2))
dom_species

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
  theme_bw() +
  theme(legend.position="bottom",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(color = guide_legend(theme = theme(legend.text = element_blank())))

dom_species / tree_height

ggsave("figure_6.png",path = here("1.LANDIS-MODEL","figures"), width = 10, height = 6.5)

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
  scale_size(range = c(0.1,2)) +
  scale_color_manual(values = onecell_colors) +
  scale_fill_manual(values = onecell_colors) +
  scale_x_continuous(breaks = c(325,375,425)) +
  scale_y_continuous(breaks = c(325,375,425)) +
  facet_wrap(~cycle, nrow=1) +
  theme_bw() +
  labs(x="X (m)",
       y = "Y (m)",
       size = "DBH (cm)") +
  theme(legend.position = "bottom",
        legend.text = element_text(face = "italic")) +
  guides(colour = guide_legend(override.aes = list(size=3), nrow=2),
         size = guide_legend(nrow=2))
stems

ggsave("figure_8.png",path = here("1.LANDIS-MODEL","figures"), width = 8.2, height = 3.5)

dom_othercell <- all_postlandis %>%
  filter(Y>=0,X>=450) %>%
  filter(Y<150,X<600)

pl_vect <- vect(dom_othercell, geom = c("X","Y"), crs = "") %>%
  rename(Species = sci_name) %>%
  mutate(cycle = factor(cycle, labels = c("10 Years\n(Spinup)",
                                          "20 Years\n(Fire cycle 1)",
                                          "30 Years\n(Fire cycle 2)",
                                          "40 Years\n(Fire cycle 3)",
                                          "50 Years\n(Fire cycle 4)")))
othercell_colors <- c("tan",colorblind_pal()(8)[c(8,3,4,2,1)])

ggplot() +
  geom_spatvector(data = pl_vect,
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
  theme_bw() +
  labs(x="X (m)",
       y = "Y (m)",
       size = "DBH (cm)") +
  theme(legend.position = "bottom",
        legend.text = element_text(face = "italic")) +
  guides(colour = guide_legend(override.aes = list(size=3), nrow=2),
         size = guide_legend(nrow=2))


ggsave("figure_8_alt.png",path = here("1.LANDIS-MODEL","figures"), width = 8.2, height = 3.5)

## figure 7
alldata_stems <- alldata_spgrp %>%
  filter(stage == "postlandis") %>%
  mutate(cycle_name = paste("Cycle", cycle)) %>%
  filter(COMMON_NAME %in% c("lodgepole pine",
                            "ponderosa pine",
                            "Douglas-fir",
                            "incense-cedar",
                            "canyon live oak")) %>%
  mutate(sci_name = factor(sci_name, levels = c("Quercus chrysolepis",
                                                "Calocedrus decurrens",
                                                "Pseudotsuga menziesii",
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
                                                "Calocedrus decurrens",
                                                "Pseudotsuga menziesii",
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
  theme_bw() +
  theme(legend.position="none")
BA_plot

ggsave("figure_7.png",path = here("1.LANDIS-MODEL","figures"), width = 8, height = 7.7)

#############
## how many trees died in each fire?
all_colors = c(colorblind_pal()(8),"tan")
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
  summarize(n_trees = n()) %>%
  pivot_wider(names_from = stage,
              values_from = n_trees) %>%
  mutate(n_dead = postlandis - postfire) %>%
  filter(fire_ID != 5) %>%
  ggplot() +
  geom_bar(aes(fire_ID, n_dead, fill = sci_name_other),
              stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(values = all_colors) +
  labs(x = "Fire cycle",
       y = "Trees killed",
       fill = "Species") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"))

ggsave("trees_killed.png",path = here("1.LANDIS-MODEL","figures"), width = 6, height = 3)

# biomass loss
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
  pivot_wider(names_from = stage,
              values_from = biomass_Mg) %>%
  mutate(biomass_loss = postlandis - postfire) %>%
  filter(fire_ID != 5) %>%
  ggplot() +
  geom_bar(aes(fire_ID, biomass_loss, fill = sci_name_other),
           stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(values = all_colors) +
  labs(x = "Fire cycle",
       y = "Biomass loss (Mg)",
       fill = "Species") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"))

ggsave("biomass_loss.png",path = here("1.LANDIS-MODEL","figures"), width = 6, height = 3)

# proportion total biomass
total_biomass <- alldata_spgrp %>%
  mutate(fire_ID = if_else(stage=="postlandis",cycle+1,cycle)) %>%
  group_by(stage, fire_ID) %>%
  summarize(total_biomass_Mg = sum(AGB_g)/1000000)

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
  left_join(total_biomass, by = join_by("fire_ID","stage")) %>%
  mutate(prop_biomass = biomass_Mg/total_biomass_Mg,
         fire_ID = fire_ID*10) %>%
  filter(stage=="postlandis") %>%
  ggplot() +
  geom_bar(aes(fire_ID, prop_biomass, fill = sci_name_other),
           stat = "identity",
           position = "stack") +
  scale_fill_manual(values = all_colors) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Years",
       y = "Percent total biomass",
       fill = "Species") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"))

ggsave("proportion_biomass.png",path = here("1.LANDIS-MODEL","figures"), width = 6, height = 3)

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

## quicfire
qf <- read.csv(here("1.LANDIS-MODEL","QF_outputs","qf_consumption.csv"))

consumption <- qf %>%
  pivot_longer(cols = 2:9,
               names_to = "var",
               values_to = "val") %>%
  mutate(strata = case_when(var %in% c("surf_cons","surf_cons_pct") ~ "Surface (<1 m)",
                            var %in% c("canopy_cons","can_cons_pct") ~ "Canopy (>1 m)",
                            var %in% c("midstory_cons","mid_cons_pct") ~ "Ladder fuels (1-7 m)",
                            var %in% c("overstory_cons","over_cons_pct") ~ "Overstory (>7 m)"),
         unit = case_when(var %in% c("surf_cons",
                                     "canopy_cons",
                                     "midstory_cons",
                                     "overstory_cons") ~ "Total (Mg)",
                          TRUE ~ "Percent")) %>%
  mutate(val = if_else(unit=="Percent",val,val/1000),
         unit = factor(unit, levels = c("Total (Mg)","Percent")),
         strata = factor(strata, levels = c("Surface (<1 m)",
                                            "Canopy (>1 m)",
                                            "Ladder fuels (1-7 m)",
                                            "Overstory (>7 m)"))) %>%
  filter(strata != "Canopy (>1 m)")

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
  theme_bw() +
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
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.minor.x = element_blank())

cons <- tot_cons / pct_cons + plot_layout(axis = "collect")
full_cons = patchwork::patchworkGrob(cons)
full_cons <- gridExtra::grid.arrange(full_cons, left = "Fine fuel consumption")

ggsave("fuel_consumption.png",
       plot = full_cons, 
       path = here("1.LANDIS-MODEL","figures"), 
       width = 5, 
       height = 3.2)
