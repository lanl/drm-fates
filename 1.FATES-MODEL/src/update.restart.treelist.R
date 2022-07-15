### Restart file outputs -------
update_restart_treelist <-
  function(sam.start,
           sam.end,
           outdir,
           VDM2FM,
           runroot,
           filebase,
           filterFile,
           finalyear,
           fire_res,
           fates_res) {
    library(ncdf4)
    library(tidyverse)

    bft <- read.table(
      file = file.path(outdir, paste0("treelist_VDM_n.plant.dat")), header = TRUE
    )
    aft <- read.table(
      file = file.path(paste0("FM2VDM/AfterFireTrees.txt")), header = FALSE
    )
    colnames(aft) <- c("fates_pft", "x", "y", "fates_height", "fates_height_to_crown_base", "fates_crown_dia",
       "height_to_widest_crown", "sizescale", "fuel_moisture_content", "bulk_density_fine_fuel", "treeid")

    # bft has grasses, but aft does not, so the difference, fire.dead.treelist, contains grasses which are also removed
    fire.dead.treelist <- bft[bft$treeid %in% setdiff(bft$treeid, aft$treeid),]
    print(paste0("No. of trees that died in fire = ", nrow(fire.dead.treelist), "; that is ", round(nrow(fire.dead.treelist)*100/length(bft$treeid), 0), "% of trees present before fire."))

    while(nrow(fire.dead.treelist) > 0) {
    # Converting each dead tree equivalent to no of trees per ha
    fire.dead.treelist$remove_nplant = 1*10000/fates_res^2      # Converting each dead tree equivalent to no of trees per ha
    fire.dead.ls <- split(fire.dead.treelist, fire.dead.treelist$nsam)
    # for each nsam in names(fire.dead.ls), go and remove n.plant in nc
    filter.arr <-
      read.table(file.path(outdir, filterFile), header = F) # on server
    sam.vec <-
      c(sam.start:sam.end)[filter.arr$V1]
    nsam <- length(sam.vec)
    # base lists
    all.sam.list <- rep(list(), nsam)
    for (i in 1:length(fire.dead.ls)) {
      sample <- sam.vec[i]
      casename <- paste0(filebase, ".", sample)
      filetag <- paste0("elm.r.", finalyear + 1, "-", "01-01-00000.nc") # a restart file at the end of final year has timestamp for the beginning of next year
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
