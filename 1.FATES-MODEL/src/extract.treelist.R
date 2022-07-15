### Restart file outputs -------
extract_treelist <- function(sam.start, sam.end, outdir, VDM2FM, runroot, filebase, var.vec.re, 
           filterFile, finalyear, fire_res, fates_res, fates_CWD_frac_twig,
           fates_c2b, leafdensity, wooddensity, sizescale_pd_df_r, grass_pft_index, 
           HYDRO, cycle_index) {
  library(ncdf4)
  library(tidyverse)
  #--------------
  # Define functions
  #--------------
  update.treeid.fun <- function(df) {
    new.df <- df[df$cycle == "new",]
    old.df <- df[df$cycle == "old",]
    # Update treeid in new.df
    # Exact matches may work for a short time-step with little growth
    exact.old.ids <- old.df$treeid[match(new.df$fates_dbh, old.df$fates_dbh)]
    new.df$treeid <- exact.old.ids  
    new.done <- new.df[!is.na(new.df$treeid),]
    
    # Left-over dataframe
    old.left <- old.df[!old.df$treeid %in% exact.old.ids,]
    new.left <- new.df[is.na(new.df$treeid),]
    
    # Because trees grow in dbh, nearest dbh, rather than exact dbh, may be the next-best predictor
    # Assuming a tree with ordered dbh as closest predictors
    # (this could be improved to find a tree with nearest dbh in the new dataset; if not already assigned a treeid.
    # Ordering both by dbh
    new.left <- new.left[order(new.left$fates_dbh),]
    old.left <- old.left[order(old.left$fates_dbh),]
    
    if (nrow(new.left) > nrow(old.left)) {
      # Trees with likely matches in the old list
      new.left$treeid[1:nrow(old.left)] <- old.left$treeid
      # New trees: Assign treeids that are not already present in the past
      remaining.rows <- c(nrow(old.left)+1):nrow(new.left)
      max.treeid <- max(old.treelist$treeid, na.rm = TRUE)
      new.df$treeid[remaining.rows] <-  max.treeid + 1:length(remaining.rows)
    } else {
      # Trees present in the old list
      new.left$treeid <- old.left$treeid[1:nrow(new.left)]
    }
    new.update.treeid <- rbind(new.done, new.left)
    return(new.update.treeid)
  }
  
  rm.duplicate.locations.fun <- function(df) {
    df$x <- runif(nrow(df), min = df$xmin[1], max = df$xmax[1]) # sample function needs size < vector of values to choose from
    df$y <- runif(nrow(df), min = df$ymin[1], max = df$ymax[1])
    # Removing grasses during duplicate location check
    df.trees <- df[df$fates_pft != as.numeric(grass_pft_index),]
    # Removing duplicates, but also reducing overlap by rounding numbers to 0.1 m (1m produces too many overlaps)
    coord <- paste0(round(df.trees$x, 1),"-", round(df.trees$y, 1))
    dupes <- which(duplicated(coord))
    if(length(dupes) > 0) {
      # Replacing duplicates with x non-overlapping with exisiting x
      nonoverlapping.set <- setdiff(seq(df.trees$xmin[1], df.trees$xmax[1], by = 0.1), round(df.trees$x, 1))
      df.trees$x[dupes] <- sample(nonoverlapping.set, length(dupes), rep = FALSE)
    }
    # Adding grasses back
    df.all <- rbind(df.trees, df[df$fates_pft == as.numeric(grass_pft_index),])
    return(df.all)
  }
  #--------------
  # End of Define functions
  #--------------
  #--------------
  # Download data
  #--------------
  
  filter.arr <-
    read.table(file.path(outdir, filterFile), header = F) # on server
  sam.vec <-
    c(sam.start:sam.end)[filter.arr$V1]
  nsam <- length(sam.vec)
  # base lists
  all.sam.list <- rep(list(), nsam) 
  for (i in 1:nsam) {
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
    nc <- nc_open(filename, write = F)
    # Write all available vars
    allvar.vec <- names(nc$var)
    write.table(
      allvar.vec,
      file = file.path(outdir, "allvar.vec.txt"),
      row.names = FALSE
    )
    res.arr <- vector("list", length = length(var.vec.re))
    for (v in 1:length(var.vec.re)) {
      var.name <- var.vec.re[v]
      names(res.arr)[v] <- var.name 
      res.arr[[v]] <- ncvar_get(nc, var.name)
    }
    nc_close(nc)
    res.all.df <- do.call(cbind.data.frame, res.arr)
    res.all.df$nsam <- i
    res.all.df$cohort.rowid <- 1:nrow(res.all.df)
    all.sam.list[[i]] <- res.all.df 
  }
  all.sam.var <- do.call(rbind, all.sam.list)
  
  #--------------
  # Make new variables
  #--------------
  all.sam.var <- all.sam.var %>%
    mutate(
      # Crown dia per plant
      fates_crown_dia = 2*(fates_cohort_area/fates_nplant/pi)^0.5,
      # Find nplants per cell, from nplants per ha
      fates_nplant.cell = fates_nplant*fates_res^2/10000, # All patch areas together amount to 1 ha or 10000 m2
      fates_height_to_crown_base = fates_height*(1-fates_crown_depth))
  
  all.sam.var <- all.sam.var %>%
    # From cohort to individual scale biomass
    mutate(fates_bagw_twig = sum(c(sapw_c_val_001, store_c_val_001, repro_c_val_001, struct_c_val_001))*fates_CWD_frac_twig*fates_c2b,
           fates_bleaf = leaf_c_val_001*fates_c2b) %>%
    select(-c(fates_crown_depth, fates_cohort_area, leaf_c_val_001, sapw_c_val_001, store_c_val_001, repro_c_val_001, struct_c_val_001))
  # Biomass weighted bulk density   
  all.sam.var <- all.sam.var %>%
    mutate(leaf_twig_bulkd = c(leafdensity*fates_bleaf + wooddensity*fates_bagw_twig)/c(fates_bleaf + fates_bagw_twig))
  
  # Biomass weighted moisture content
  if (HYDRO == TRUE) {
    moisture <- read.table(file = file.path(VDM2FM, "livefuel.moisture.txt"), header = TRUE)
    all.sam.var <- all.sam.var %>%
      left_join(moisture, by = c("nsam", "fates_pft")) %>%
      mutate(fuel_moisture_content = c(FATES_LTH_SCPF*fates_bleaf + FATES_STH_SCPF*fates_bagw_twig)/c(fates_bleaf + fates_bagw_twig)) %>%
      select(-FATES_LTH_SCPF, -FATES_STH_SCPF)
  } else {
    all.sam.var <- all.sam.var %>%
      mutate(fuel_moisture_content = 0.5)
  }
  
  #--------------
  ## From nplants per size-pft to individual plants
  ##--------------
  trees.whole <- all.sam.var %>%
    # Selecting only whole trees to present to the Fire model
    # This could be made complex by grouping trees into size class bins and checkign if the size class makes it to a whole tree (nplant.400 = 1)
    mutate(fates_nplant.cell = as.integer(round(fates_nplant.cell, 0))) %>%
    subset(fates_nplant.cell >= 1) %>%
    # repeat each plant number of plants
    uncount(fates_nplant.cell)
  trees.whole$treeid <- 1:nrow(trees.whole)
  
  #--------------
  ## Assign x-y location for each FATES simulation (aka patch). Assume simulation index, nsam, increases from East to West, then increases northwards.
  #--------------
  
  nsam_side = sam.end^0.5
  cell.xy <- data.frame(nsam = c(1:sam.end),
                        xmin = rep(fates_res*c(1:nsam_side-1) + 1, times = nsam_side),
                        ymin = rep(fates_res*c(1:nsam_side-1) + 1, each = nsam_side)) %>%
    mutate(xmax = xmin - 1 + fates_res, ymax = ymin - 1 + fates_res)
  
  trees.whole <- trees.whole %>%
    left_join(cell.xy, by = "nsam") 
  
  #------------------------------------
  # Assign tree locations for new trees
  #------------------------------------
  
  # Get tree locations from last FATES run
  if (cycle_index > 0) {
    old.treelist <- read.table(
      file = file.path(VDM2FM, paste0("treelist_VDM_n.plant.", cycle_index - 1, ".dat")), header = TRUE)

    new.id.df <- trees.whole[, c("nsam", "treeid", "fates_pft", "fates_dbh")]
    old.id.df <- old.treelist[, c("nsam", "treeid", "fates_pft", "fates_dbh")]
    new.id.df$cycle <- "new"; old.id.df$cycle <- "old"
    id.df <- rbind(new.id.df, old.id.df)
    new.id.df$rownum <- 1:nrow(new.id.df)
    
    id.df.ls <- split(id.df, f=list(id.df$nsam, id.df$fates_pft), drop = TRUE)
    
    # Within each nsam by pft dataframe, update treeid for new treelist
    
    update.treeid.ls <- lapply(id.df.ls, update.treeid.fun)
    
    update.treeid.df <- dplyr::bind_rows(update.treeid.ls)
    # because new treeids were created within nsam loop, there may be duplicates
    trees.whole$treeid <- update.treeid.df$treeid
  }
  
  # ----
  ## Removing coordinate duplicates within each FATES sim, since min max differ for each sim
  # ---
  
  treelist.ls <- split(trees.whole, f=list(trees.whole$nsam), drop = TRUE)
  treelist.ls <- lapply(treelist.ls, rm.duplicate.locations.fun)
  # Checking if any tree location duplicates remain
  treelist <- dplyr::bind_rows(treelist.ls)
  treelist.no.grass <- treelist[treelist$fates_pft != as.numeric(grass_pft_index),]
  coord <- paste0(round(treelist.no.grass$x, 1),"-", round(treelist.no.grass$y, 1))
  dupes <- which(duplicated(coord))
  print(paste0("Remaining coordinate duplicates within 1m distance = ", length(dupes)))
  treelist <- treelist %>%
    mutate(height_to_widest_crown = fates_height_to_crown_base, # not ideal for non-pines
           bulk_density_fine_fuel = leaf_twig_bulkd) %>%
    left_join(sizescale_pd_df_r, by = "fates_pft") %>%
    select(c(fates_pft, x, y, fates_height, fates_height_to_crown_base, fates_crown_dia,
             height_to_widest_crown, sizescale, fuel_moisture_content,
             bulk_density_fine_fuel, treeid, nsam, fates_nplant, cohort.rowid, fates_dbh))
  treelist <- as.data.frame(sapply(treelist, as.numeric))
  treelist.no.grass <- treelist[which(treelist$fates_pft != as.numeric(grass_pft_index)),]
  treelist.VDM2FM <- subset(treelist.no.grass, select = -c(nsam, fates_nplant, cohort.rowid, fates_dbh))

  write.table(
    treelist,
    file = file.path(VDM2FM, paste0("treelist_VDM_n.plant.", cycle_index, ".dat")),
    row.names = FALSE
  )
  
  if(!dir.exists(VDM2FM)) {dir.create(VDM2FM)}
  write.table(
    treelist.VDM2FM,
    file = file.path(VDM2FM, paste0("treelist_VDM.dat")),
    row.names = FALSE, col.names = FALSE
  )
  write.table(
    treelist.VDM2FM,
    file = file.path(VDM2FM, paste0("treelist_VDM.", cycle_index, ".dat")),
    row.names = FALSE, col.names = FALSE
  )
  return(TRUE)
}


