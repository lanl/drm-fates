library(here)
library(tidyverse)
library(terra)
library(ggthemes)
library(tidyterra)
library(collapse)
library(patchwork)
library(cowplot)

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

alldata_count <-  alldata_spgrp %>%
  group_by(cycle,stage,MapCode) %>%
  summarize(avg_height = mean(HT_m),
            func_grp = fmode(func_grp, na.rm=T),
            dom_count = fmode(sci_name,na.rm=T))

alldata_biomass <- alldata_spgrp %>%
  group_by(cycle,stage,MapCode,sci_name) %>%
  summarize(biomass_sum = sum(AGB_g, na.rm=T)) %>%
  ungroup() %>%
  group_by(cycle,stage,MapCode) %>%
  slice_max(biomass_sum, n=1, with_ties=F) %>%
  rename(dom_biomass=sci_name) %>%
  select(cycle,MapCode,dom_biomass)

alldata_BA <- alldata_spgrp %>%
  group_by(cycle,stage,MapCode,sci_name) %>%
  summarize(DBH_sum = sum(DBH, na.rm=T)) %>%
  ungroup() %>%
  group_by(cycle,stage,MapCode) %>%
  slice_max(DBH_sum, n=1, with_ties=F) %>%
  rename(dom_BA=sci_name) %>%
  select(cycle,MapCode,dom_BA)

alldata_dom <- left_join(alldata_count,alldata_biomass, by = join_by("cycle","stage","MapCode"))
alldata_dom <- left_join(alldata_dom,alldata_BA, by=join_by("cycle","stage","MapCode"))

alldata_fact <- alldata_dom %>%
  mutate(func_grp = as.integer(factor(func_grp)),
         dom_count = as.integer(factor(dom_count)),
         dom_biomass = as.integer(factor(dom_biomass)),
         dom_BA = as.integer(factor(dom_BA)))

t0_pl <- cycle_rst(0, "postlandis", OC_df, alldata_fact, alldata_dom)
t1_pl <- cycle_rst(1, "postlandis", OC_df, alldata_fact, alldata_dom)
t2_pl <- cycle_rst(2, "postlandis", OC_df, alldata_fact, alldata_dom)
t3_pl <- cycle_rst(3, "postlandis", OC_df, alldata_fact, alldata_dom)
t4_pl <- cycle_rst(4, "postlandis", OC_df, alldata_fact, alldata_dom)

###################
treelists_vect <- vect(alldata_spgrp, geom=c("X","Y"))

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
names(BA_rasters) <- c("Cycle 0", "Cycle 1","Cycle 2","Cycle 3","Cycle 4")
###################

# Recreate top row of zach's plot
basal_area <- c(
  t0_pl$dom_BA,
  t1_pl$dom_BA,
  t2_pl$dom_BA,
  t3_pl$dom_BA,
  t4_pl$dom_BA
)
names(basal_area) <- c("10 Years\n(Spinup)",
                       "20 Years\n(Fire cycle 1)",
                       "30 Years\n(Fire cycle 2)",
                       "40 Years\n(Fire cycle 3)",
                       "50 Years\n(Fire cycle 4)")
basal_area <- basal_area %>%
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

all_colors <- c(
  "#B0E0E6",
  "#FFA07A",
  "#FF7F50",
  "#FF6347",
  "#FF4500",
  "#87CEEB",
  "#DC143C",
  "#B22222",
  "#4682B4",
  "#8B0000",
  "#00008B"
)


perim1 <- read_perimeters(1)
perim2 <- read_perimeters(2)
perim3 <- read_perimeters(3)
perim4 <- read_perimeters(4)
perimeters <- rbind(perim1,perim2,perim3,perim4)
perimeters <- perimeters %>%
  mutate(lyr = case_when(cycle==1 ~ "20 Years\n(Fire cycle 1)",
                           cycle==2 ~ "30 Years\n(Fire cycle 2)",
                           cycle==3 ~ "40 Years\n(Fire cycle 3)",
                           TRUE ~ "50 Years\n(Fire cycle 4)"))

dom_species <- ggplot() +
  geom_spatraster(data=basal_area) +
  geom_spatvector(data=perimeters, color = "red", fill="black", alpha = 0.2) +
  scale_fill_colorblind() +
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
height <- c(
  t0_pl$avg_height,
  t1_pl$avg_height,
  t2_pl$avg_height,
  t3_pl$avg_height,
  t4_pl$avg_height
)
names(height) <- c("10 Years\n(Spinup)",
                   "20 Years\n(Fire cycle 1)",
                   "30 Years\n(Fire cycle 2)",
                   "40 Years\n(Fire cycle 3)",
                   "50 Years\n(Fire cycle 4)")
ext(height) = ext(c(0,750,0,750))

tree_height <- ggplot() +
  geom_spatraster(data=height) +
  geom_spatvector(data=perimeters, color = "red", fill=NA) +
  scale_fill_distiller(palette = 2, direction=1, na.value="white") +
  facet_wrap(~lyr, nrow=1) +
  labs(x = "X (m)",
       y = "Y (m)",
       fill = "Average Tree Height (m)") +
  theme_bw() +
  theme(legend.position="bottom",
        strip.background = element_blank(),
        strip.text = element_blank())

dom_species / tree_height

ggsave("figure_6.png",path = here("1.LANDIS-MODEL","figures"), width = 8, height = 5.8)

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
onecell_colors <- colorblind_pal()(8)[c(2,3,4)]

stems <- ggplot() +
  geom_spatvector(data = pl_vect,
                aes(color = Species, 
                    fill = Species,
                    size = DBH),
                shape = 21,
                alpha = 0.5) +
  scale_size(range = c(0.5,4)) +
  # scale_color_manual(values = onecell_colors) +
  # scale_fill_manual(values = onecell_colors) +
  scale_x_continuous(breaks = c(325,375,425)) +
  scale_y_continuous(breaks = c(325,375,425)) +
  facet_wrap(~cycle, nrow=1) +
  theme_bw() +
  labs(x="X (m)",
       y = "Y (m)",
       size = "DBH (cm)") +
  theme(legend.position = "bottom",
        legend.text = element_text(face = "italic"))
stems

ggsave("figure_8.png",path = here("1.LANDIS-MODEL","figures"), width = 8.2, height = 3)

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



