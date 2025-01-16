library(here)
library(tidyverse)
library(terra)
library(ggthemes)
library(tidyterra)
library(collapse)
library(patchwork)

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

# Recreate top row of zach's plot
basal_area <- c(
  t0_pl$dom_BA,
  t1_pl$dom_BA,
  t2_pl$dom_BA,
  t3_pl$dom_BA,
  t4_pl$dom_BA
)
names(basal_area) <- c("10 Years","20 Years","30 Years","40 Years","50 Years")
basal_area <- basal_area %>%
  mutate(`10 Years` = factor(`10 Years`, 
                             levels = c("Quercus chrysolepis",
                                        "Pseudotsuga menziesii",
                                        "Calocedrus decurrens",
                                        "Pinus ponderosa",
                                        "Pinus contorta",
                                        "Taxus brevifolia")),
         `20 Years` = factor(`20 Years`, 
                             levels = c("Quercus chrysolepis",
                                        "Pseudotsuga menziesii",
                                        "Calocedrus decurrens",
                                        "Pinus ponderosa",
                                        "Pinus contorta",
                                        "Taxus brevifolia")),
         `30 Years` = factor(`30 Years`, 
                             levels = c("Quercus chrysolepis",
                                        "Pseudotsuga menziesii",
                                        "Calocedrus decurrens",
                                        "Pinus ponderosa",
                                        "Pinus contorta",
                                        "Taxus brevifolia")),
         `40 Years` = factor(`40 Years`, 
                             levels = c("Quercus chrysolepis",
                                        "Pseudotsuga menziesii",
                                        "Calocedrus decurrens",
                                        "Pinus ponderosa",
                                        "Pinus contorta",
                                        "Taxus brevifolia")),
         `50 Years` = factor(`50 Years`, 
                             levels = c("Quercus chrysolepis",
                                        "Pseudotsuga menziesii",
                                        "Calocedrus decurrens",
                                        "Pinus ponderosa",
                                        "Pinus contorta",
                                        "Taxus brevifolia")))

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


dom_species <- ggplot() +
  geom_spatraster(data=basal_area) +
  scale_fill_colorblind() +
  facet_wrap(~lyr, nrow=1) +
  coord_equal() +
  theme_bw() +
  theme(legend.position="bottom",
        legend.text = element_text(face="italic"),
        legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
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
names(height) <- c("10 Years","20 Years","30 Years","40 Years","50 Years")

tree_height <- ggplot() +
  geom_spatraster(data=height) +
  scale_fill_distiller(palette = 2, direction=1, na.value="white") +
  facet_wrap(~lyr, nrow=1) +
  coord_equal() +
  labs(fill = "Average Tree Height (m)") +
  theme_bw() +
  theme(legend.position="bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank())

dom_species / tree_height

ggsave("figure_6_ggplot.png",path = here("LANDIS_DRM_maps"), width = 8, height = 7.7)


# figure 7
overlap_stems <- overlap_spgrp %>%
  mutate(cycle_name = paste("Cycle", cycle)) %>%
  filter(COMMON_NAME %in% c("canyon live oak",
                            "Douglas-fir",
                            "ponderosa pine",
                            "lodgepole pine")) %>%
  ggplot() +
  geom_histogram(aes(DBH), bins = 12) +
  facet_grid(sci_name ~ cycle_name, scales = "free_y") +
  theme_bw() +
  theme()
overlap_stems
