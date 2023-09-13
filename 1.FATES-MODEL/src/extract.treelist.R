
# Notes
# 1. Refer to "DRM 03 (Feb04, 2023).docx" for details

### restart file outputs -------
extract_treelist <- function(sam.start, sam.end, outdir, VDM2FM, runroot, filebase, var.vec.re, 
           filterFile, finalyear, fire_res, fates_res, fates_CWD_frac_twig,
           fates_c2b, leafdensity, wooddensity, sizescale_pd_df_r, grass_pft_index, 
           HYDRO, cycle_index) {

  #CSXM: load packages
  library(ncdf4)
  library(tidyverse)

  # print(paste("in R, grass_pft_index", grass_pft_index)) #ASXM

  #----------------
  # Functions (BGN)
  #----------------
  update.treeid.fun <- function(df) {
    new.df <- df[df$cycle == "new",]
    old.df <- df[df$cycle == "old",]
    # Update treeid in new.df
    # Exact matches may work for a short time-step with little growth
    exact.old.ids <- old.df$treeid[match(new.df$fates_dbh, old.df$fates_dbh)]
    new.df$treeid <- exact.old.ids  
    new.done <- new.df[!is.na(new.df$treeid),]
    if(length(which(is.na(new.df$treeid))) > 0) {
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
    } else {
      new.update.treeid <- new.done
    }
    return(new.update.treeid)
  }
  
  # commented out by SXM (BGN)
  # rm.duplicate.locations.fun <- function(df) {
  #   df$x <- runif(nrow(df), min = df$xmin[1], max = df$xmax[1]) # sample function needs size < vector of values to choose from
  #   df$y <- runif(nrow(df), min = df$ymin[1], max = df$ymax[1])
  #   # Removing grasses during duplicate location check
  #   df.trees <- df[df$fates_pft != as.numeric(grass_pft_index),]
  #   # Removing duplicates, but also reducing overlap by rounding numbers to 0.1 m (1m produces too many overlaps)
  #   coord <- paste0(round(df.trees$x, 1),"-", round(df.trees$y, 1))
  #   dupes <- which(duplicated(coord))
  #   if(length(dupes) > 0) {
  #     # Replacing duplicates with x non-overlapping with exisiting x
  #     nonoverlapping.set <- setdiff(seq(df.trees$xmin[1], df.trees$xmax[1], by = 0.1), round(df.trees$x, 1))
  #     df.trees$x[dupes] <- sample(nonoverlapping.set, length(dupes), rep = FALSE)
  #   }
  #   # Adding grasses back
  #   df.all <- rbind(df.trees, df[df$fates_pft == as.numeric(grass_pft_index),])
  #   return(df.all)
  # }
  # commented out by SXM (END)

  # Added by SXM (BGN)
  # Refn: /usr/projects/climate/xiaoming/Analysis/zTech/scrpts.R/Scrpts.Algorithm/algorithm_2Dxy_SXM.R
  rm.duplicate.locations.fun <- function(df) {
    dBug = TRUE
    df$x <- runif(nrow(df), min = df$xmin[1], max = df$xmax[1]) 
    df$y <- runif(nrow(df), min = df$ymin[1], max = df$ymax[1])
    if(dBug) {
      cat("\n\nrm.duplicate.locations.fun: dBug01 \n")
      print(paste("df$xmin[1]=", df$xmin[1], "df$xmax[1]=", df$xmax[1]))
      print(paste("nrow(df)=", nrow(df)))
    }

    # remove grasses
    df.trees <- df[!(df$fates_pft %in% grass_pft_index), ]
    if(dBug) {
      cat("\n\nrm.duplicate.locations.fun: dBug02 \n")
      print(paste("nrow(df.trees)=", nrow(df.trees)))
    }

    # uniformly distribute df.trees in the spatial domain
    # step01: setups for random num generation to make results reproducible (ensemble effects has already been considered in simus)
    seedNum = 123
 
    # step02: refine iteratively to avoid overlapped locations (locations wt. a distance to other locations less than 0.1 m is defined as overlapped)
    coord <- paste0(round(df.trees$x, 1), "-", round(df.trees$y, 1))
    df.trees.dupe = df.trees[duplicated(coord),  ]
    df.trees      = df.trees[!duplicated(coord), ]
    while(nrow(df.trees.dupe) > 0) {
      # regenerate nrow(df.trees.dupe) num of uniformaly distributed (x,y) locations to replace the duplicated (x,y) locations in df.trees
      set.seed(seedNum)
      df.trees.dupe <- df.trees.dupe %>%
        mutate(
          x = runif(nrow(df.trees.dupe), df$xmin[1], df$xmax[1]), 
          y = runif(nrow(df.trees.dupe), df$ymin[1], df$ymax[1]))
      df.trees = rbind(df.trees, df.trees.dupe)
      if(dBug) {
        cat("\n\nrm.duplicate.locations.fun: dBug03 within the while-loop \n")
        print(paste("nrow(df.trees.dupe) = ", nrow(df.trees.dupe)))
        print(paste("nrow(df.trees)      = ", nrow(df.trees))     )
      }

      # identify duplicated (x,y) and remove them
      coord <- paste0(round(df.trees$x, 1), "-", round(df.trees$y, 1))
      df.trees.dupe = df.trees[duplicated(coord), ]
      if(nrow(df.trees) > 500000 && nrow(df.trees.dupe) <= 5) { # allow up to 5 trees to overlap wt. other trees if totTrees is greater than 500000
        break
      }
      df.trees = df.trees[!duplicated(coord), ]
    }

    # add grasses back
    df.all <- rbind(df.trees, df[df$fates_pft %in% grass_pft_index,])

    return(df.all)
  }
  # Added by SXM (END)
  #----------------
  # Functions (END)
  #----------------

  #+++++++++++++++++++
  # Read FATES Output 
  #+++++++++++++++++++
  filter.arr <-
    read.table(file.path(outdir, filterFile), header=F) # on server #CSXM: list
  sam.vec <-
    c(sam.start:sam.end)[filter.arr$V1] #CSXM: only include completed FATES simus; filter.arr$V1 = [TRUE, TRUE, FALSE, TRUE] indicate the 3rd FATES simu is not completed and not inlcuded
  nsam <- length(sam.vec)

  # base lists
  all.sam.list <- rep(list(), nsam) #CSXM: list is a relatively complex structure in R, dfrnt from Python 
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
    nc <- nc_open(filename, write=F)

    # write all available vars in FATES restarting files into allvar.vec.txt
    allvar.vec <- names(nc$var)
    write.table(
      allvar.vec,
      file = file.path(outdir, "allvar.vec.txt"),
      row.names = FALSE
    )
    res.arr <- vector("list", length = length(var.vec.re)) #CSXM: a vector of lists
    for (v in 1:length(var.vec.re)) {
      var.name <- var.vec.re[v]
      names(res.arr)[v] <- var.name 
      res.arr[[v]] <- ncvar_get(nc, var.name)
    }
    nc_close(nc)
    res.all.df <- do.call(cbind.data.frame, res.arr) #CSXM: combine vars in res.arr as a data frame and store it in res.all.df
    res.all.df$nsam <- i                             #CSXM: in each row of res.all.df, add a column named as nsam
    res.all.df$cohort.rowid <- 1:nrow(res.all.df)    #CSXM: in each row of res.all.df, add a column named as cohort.rowid
    all.sam.list[[i]] <- res.all.df                  #CSXM: store res.all.df as the ith element in the all.sam.list (a list)
  }
  all.sam.var <- do.call(rbind, all.sam.list)

  #++++++++++++++++++++
  # Make New Variables
  #++++++++++++++++++++
  all.sam.var <- all.sam.var %>% #CSXM: the pipe operator that fowards a value, or the result of an experssion, into the next function call/expression (i.e., pass the left hand side of the operator to the FIRST argument of the right hand side of the operator)
    mutate(
      # crown dia per plant
      fates_crown_dia = 2*(fates_cohort_area/fates_nplant/pi)^0.5,

      # find nplants per cell, from nplants per ha
      fates_nplant.cell = fates_nplant*fates_res^2/10000, # all patch areas together amount to 1 ha or 10000 m2

      #CSXM: distance from tree base to crown base
      # fates_height_to_crown_base = fates_height*(1-fates_crown_depth)) # commented out by SXM
      fates_height_to_crown_base = fates_height*(1.0-fates_crown_depth/100.0)) # ASXM: fates_crown_depth in the restart-.r. files is called fraction, but actually in %

  # from cohort to individual scale biomass 
  all.sam.var <- all.sam.var %>%
    mutate(fates_bagw_twig = sum(c(sapw_c_val_001, store_c_val_001, repro_c_val_001, struct_c_val_001))*fates_CWD_frac_twig*fates_c2b, #QSXM01: bagw means biomass of above ground wood twig and where is this formula from?
           fates_bleaf = leaf_c_val_001*fates_c2b) %>%
    select(-c(fates_crown_depth, fates_cohort_area, leaf_c_val_001, sapw_c_val_001, store_c_val_001, repro_c_val_001, struct_c_val_001))

  # biomass weighted bulk density 
  all.sam.var <- all.sam.var %>%
    mutate(leaf_twig_bulkd = c(leafdensity*fates_bleaf + wooddensity*fates_bagw_twig)/c(fates_bleaf + fates_bagw_twig)) #CSXM: weighted based on leaf and above gound wood twig

  # biomass weighted moisture content #CSXM: HYDRO is never TRUE in DRM, such that its fuel_moisture_content = 0.5 always
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
  
  #++++++++++++++++++++++++++++++++++++++++++++++++
  # From nplants per size-pft to Individual Plants
  #++++++++++++++++++++++++++++++++++++++++++++++++
  #CSXM: in trees.whole, each line is a descripiton of each tree based on FATES output (the so-called individual plants)
  trees.whole <- all.sam.var %>%
    # selecting only whole trees to present to the fire model
    # this could be made complex by grouping trees into size class bins and checking if the size class makes it to a whole tree (nplant.400 = 1)
    mutate(fates_nplant.cell = as.integer(round(fates_nplant.cell, 0))) %>%
    subset(fates_nplant.cell >= 1) %>%
    # repeat each plant number of plants #CSXM: expand the data frame by repeating each plant-line num-of-plants times (num-of-plants have been rounded as integer two lines above)
    uncount(fates_nplant.cell)
  trees.whole$treeid <- 1:nrow(trees.whole)

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Assign x-y location for each FATES simulation (aka patch) 
  # (assume simulation index, nsam, increases from East to West, then increases northward)
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #CSXM: indexing the domain (assumed 1 meter btw each indexed point in the x and y directions)
  # cell.xy:
  #   nsam xmin ymin xmax ymax
  # 1    1    1    1  200  200
  # 2    2  201    1  400  200
  # 3    3    1  201  200  400
  # 4    4  201  201  400  400
  nsam_side = sam.end^0.5 #CSXM: the number of point runs must be an integer after the square root operation
  cell.xy <- data.frame(nsam = c(1:sam.end),
                        xmin = rep(fates_res*c(1:nsam_side-1) + 1, times = nsam_side),
                        ymin = rep(fates_res*c(1:nsam_side-1) + 1, each = nsam_side)) %>%
    mutate(xmax = xmin - 1 + fates_res, ymax = ymin - 1 + fates_res)

  #CSXM: each invidual plant is tied to their simulation domain (take each point run as a domain with one grid cell)
  trees.whole <- trees.whole %>%
    left_join(cell.xy, by="nsam") 

  #+++++++++++++++++++++++++++++++++++++
  # Assign tree locations for new trees
  #+++++++++++++++++++++++++++++++++++++
# During the second time of reading this script, we stopped here on Apr. 5 2023 
# because do not think further information here is needed for us to figure out the PDF of tree height that is to be used in FATES calibration


  # get tree locations from last FATES run (#CSXM: last means the most recent previous)
  if (cycle_index > 0) {
    old.treelist <- read.table(
      file = file.path(VDM2FM, paste0("treelist_VDM_n.plant.", cycle_index-1, ".dat")), header = TRUE)

    new.id.df <- trees.whole[, c("nsam", "treeid", "fates_pft", "fates_dbh")]
    old.id.df <- old.treelist[, c("nsam", "treeid", "fates_pft", "fates_dbh")]
    new.id.df$cycle <- "new"; old.id.df$cycle <- "old"
    id.df <- rbind(new.id.df, old.id.df)
    id.df.ls <- split(id.df, f=list(id.df$nsam, id.df$fates_pft), drop = TRUE)
    # within each nsam by pft dataframe, update treeid for new treelist
    
    update.treeid.ls <- lapply(id.df.ls, update.treeid.fun)
    update.treeid.df <- dplyr::bind_rows(update.treeid.ls)
    # because new treeids were created within nsam loop, there may be duplicates
    trees.whole$treeid <- update.treeid.df$treeid
  }

  #------------------------------------------------------------------------------------------
  # Removing coordinate duplicates within each FATES simu, since min max differ for each simu
  #------------------------------------------------------------------------------------------
  # group according to simus (sam.start to sam.end)
  cat("\n\nBEFFORE grouping AND BEFORE getting into rm.duplicate.locations.fun \n") # ASXM
  print(paste("nrow(trees.whole)", nrow(trees.whole))) # ASXM
  treelist.ls <- split(trees.whole, f=list(trees.whole$nsam), drop = TRUE) # treelist.ls may be the same with f=list(trees.whole$nsam) or f=trees.whole$nsam (to check)

  # distribute trees uniformly in the domain (keep using rm.duplicate.locations.fun as the function name, although it is confusing)
  # CSXM: because treelist.ls is grouped for simus (four in total), the following wil call the rm.duplicate.locations four timess
  treelist.ls <- lapply(treelist.ls, rm.duplicate.locations.fun)

  #*********************************************************************************************************
  # CSXM: Roughly stopped reading here during the last visit in Feb 2023
  # Reading notes: /usr/projects/climate/xiaoming/Analysis/zTech/scrpts.R/Scrpts.DRM/drm_extract_treelist.R
  #*********************************************************************************************************

  # Checking if any tree location duplicates remain
  treelist <- dplyr::bind_rows(treelist.ls)
  #CSXM: for multiPFTs, wrt as "treelist.no.grass <- treelist[!(treelist$fates_pft  grass_pft_index),]"
  treelist.no.grass <- treelist[treelist$fates_pft != as.numeric(grass_pft_index),]
  coord <- paste0(round(treelist.no.grass$x, 1),"-", round(treelist.no.grass$y, 1))
  dupes <- which(duplicated(coord))
  # print(paste0("Remaining coordinate duplicates within 1m distance = ", length(dupes))) # commented out by SXM
  print(paste0("Remaining coordinate duplicates within 0.1 m distance = ", length(dupes)))

  # CSXM: in the following from select(c(fates_pft, ..., fates_dbh)), the first 10 fields are in the same order as in update.restart.treelist.R, although not sure it is actually needed
  treelist <- treelist %>%
    mutate(height_to_widest_crown = fates_height_to_crown_base, # not ideal for non-pines
           bulk_density_fine_fuel = leaf_twig_bulkd) %>% 
    left_join(sizescale_pd_df_r, by = "fates_pft") %>%
    select(c(fates_pft, x, y, fates_height, fates_height_to_crown_base, fates_crown_dia,   # CSXM: column 1-6 (fates_height_to_crown_base is actually the bottom of canopy)
             height_to_widest_crown, sizescale, fuel_moisture_content,                     # CSXM: column 7-9
             bulk_density_fine_fuel, treeid, nsam, fates_nplant, cohort.rowid, fates_dbh)) # CSXm: column 10-15
  treelist <- as.data.frame(sapply(treelist, as.numeric))
  #CSXM: for multiPFTs, wrt as "treelist.no.grass <- treelist[which(!(treelist$fates_pft  grass_pft_index)),]"
  treelist.no.grass <- treelist[which(treelist$fates_pft != as.numeric(grass_pft_index)),] 
  # treelist.VDM2FM <- subset(treelist.no.grass, select = -c(nsam, fates_nplant, cohort.rowid, fates_dbh)) # commented out by SXM
  treelist.VDM2FM <- subset(treelist.no.grass, select = -c(nsam, fates_nplant, cohort.rowid, fates_dbh))

  write.table(
    treelist,
    file = file.path(VDM2FM, paste0("treelist_VDM_n.plant.dat")),
    row.names = FALSE
  )

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
