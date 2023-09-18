### Restart file outputs -------
update_restart_treelist <-
  function(sam.start,
           sam.end,
           outdir,
           VDM2FM,
           FM2VDM,
           runroot,
           filebase,
           filterFile,
           finalyear,
           fire_res,
           fates_res,
	   cycle_index,
           grass_pft_index) {
    library(ncdf4)
    library(tidyverse)

    bft <- read.table(
      # file = file.path(paste0("VDM2FM/treelist_VDM_n.plant.dat")), header = TRUE # commented out by SXM
      file = file.path(paste0(VDM2FM, "/treelist_VDM_n.plant.dat")), header = TRUE # ASXM
    )
    aft <- read.table(
      # file = file.path(paste0("FM2VDM/AfterFireTrees.txt")), header = FALSE # commented out by SXM
      file = file.path(paste0(FM2VDM, "/AfterFireTrees.txt")), header = FALSE
    )
    colnames(aft) <- c("fates_pft", "x", "y", "fates_height", "fates_height_to_crown_base", "fates_crown_dia",
        "height_to_widest_crown", "sizescale", "fuel_moisture_content", "bulk_density_fine_fuel", "treeid")

    # bft has grasses, but aft does not, so the difference, fire.dead.plantlist, contains grasses which are also removed
    fire.dead.plantlist <- bft[bft$treeid %in% setdiff(bft$treeid, aft$treeid),]

    print(paste0("No. of plants that died in fire = ", nrow(fire.dead.plantlist), "; that is ", round(nrow(fire.dead.plantlist)*100/length(bft$treeid), 0), "% of plants present before fire."))

    #CSXM: for multiPFTs, wrt as "bft.trees <- bft[!(bft$fates_pft %in% grass_pft_index), ]"
    bft.trees <- bft[bft$fates_pft != as.numeric(grass_pft_index),] 

    #CSXM: for multiPFTs, wrt as "fire.dead.treelist <- fire.dead.plantlist[!(fire.dead.plantlist$fates_pft %in% grass_pft_index), ]"
    fire.dead.treelist <- fire.dead.plantlist[fire.dead.plantlist$fates_pft != as.numeric(grass_pft_index),] 

    print(paste0("No. of trees that died in fire = ", nrow(fire.dead.treelist), "; that is ", round(nrow(fire.dead.treelist)*100/length(bft.trees$treeid), 0), "% of trees present before fire."))

    while(nrow(fire.dead.plantlist) > 0) {
    # Converting each dead tree equivalent to no of trees per ha
    fire.dead.plantlist$remove_nplant = 1*10000/fates_res^2      # Converting each dead tree equivalent to no of trees per ha
    fire.dead.ls <- split(fire.dead.plantlist, fire.dead.plantlist$nsam)
    # for each nsam in names(fire.dead.ls), go and remove n.plant in nc
    filter.arr <-
      # read.table(file.path(outdir, filterFile), header = F) # on server # commented out by SXM
      read.table(file.path(paste0(outdir, "/", filterFile)), header = F) # ASXM
    sam.vec <-
      c(sam.start:sam.end)[filter.arr$V1]
    nsam <- length(sam.vec)

    # base lists
    all.sam.list <- rep(list(), nsam)
    for (i in 1:length(fire.dead.ls)) {
      sample <- sam.vec[i]
      casename <- paste0(filebase, ".", sample)
      if (cycle_index==0) {
	filetag <- paste0("elm.r.", finalyear, "-", "01-01-01800.nc") # at the end of the first month
      } else {
	filetag <- paste0("elm.r.", finalyear + 1, "-", "01-01-01800.nc") # a restart file at the end of final year has timestamp for the beginning of next year
      }
      filename <-
        paste0(runroot,
               "/",
               casename,
               "/run/",
               casename,
               ".",
               filetag)
      nc <- nc_open(filename, write = TRUE)
      old.n.plant <- ncvar_get(nc, "fates_nplant")
      nsam.df <- fire.dead.ls[[names(fire.dead.ls)[i]]]
      post_fire_n.plant <- old.n.plant
      post_fire_n.plant[nsam.df$cohort.rowid] <- old.n.plant[nsam.df$cohort.rowid] - nsam.df$remove_nplant
      ncdf4::ncvar_put(nc, "fates_nplant", post_fire_n.plant)
      nc_close(nc)
    }
    # Test substitution for the last file:
    nc <- nc_open(filename, write = F)
    ## value in the file
    file.value <- ncvar_get(nc, "fates_nplant")
    ## value that was to be substituted
    desired.value <- post_fire_n.plant[nsam.df$cohort.rowid]
    ncdf4::nc_close(nc)
    ## make sure the two do match
    match <- all(file.value[nsam.df$cohort.rowid] == desired.value)

    return(match)
    }
  }
