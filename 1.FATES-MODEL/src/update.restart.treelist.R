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
    #library(tidyverse)

    # ******************
    # Find trees and grasses that died. 
    # Each has a corresponding cohort ID they come from via Columnname cohort.rowid
    # ******************

    bft <- read.table(
      file = file.path(paste0(VDM2FM, "/treelist_VDM_n.plant.", cycle_index, ".dat")), header = TRUE # ASXM
    )
    aft <- read.table(
      file = file.path(paste0(FM2VDM, "/AfterFireTrees.", cycle_index + 1, ".txt")), header = FALSE
    )
    colnames(aft) <- c("fates_pft", "x", "y", "fates_height", "fates_height_to_crown_base", "fates_crown_dia",
        "height_to_widest_crown", "sizescale", "fuel_moisture_content", "bulk_density_fine_fuel", "treeid")
    # list all the plants that died by taking the difference between before fire plant list (includes grasses) and after fire tree list.
    fire.dead.plantlist <- bft[bft$treeid %in% setdiff(bft$treeid, aft$treeid),] # so fire.dead.plantlist contains grasses, because bft has grasses, but aft does not 

    print(paste0("No. of plants that died in fire = ", nrow(fire.dead.plantlist), "; that is ", round(nrow(fire.dead.plantlist)*100/length(bft$treeid), 0), "% of plants present before fire."))

    #CSXM: for multiPFTs, wrt as "bft.trees <- bft[!(bft$fates_pft %in% grass_pft_index), ]"
    bft.trees <- bft[bft$fates_pft != as.numeric(grass_pft_index),] 

    #CSXM: for multiPFTs, wrt as "fire.dead.treelist <- fire.dead.plantlist[!(fire.dead.plantlist$fates_pft %in% grass_pft_index), ]"
    fire.dead.treelist <- fire.dead.plantlist[fire.dead.plantlist$fates_pft != as.numeric(grass_pft_index),] 

    print(paste0("No. of trees that died in fire = ", nrow(fire.dead.treelist), "; that is ", round(nrow(fire.dead.treelist)*100/length(bft.trees$treeid), 0), "% of trees present before fire."))

    # ******************
    # Convert dead trees and grasses into cohort densities
    # Sum dead plant densities per cohort ID
    # ******************

    while(nrow(fire.dead.plantlist) > 0) {
    # Converting each dead tree equivalent to no of trees per ha
    fire.dead.plantlist$remove_nplant = 1*10000/fates_res^2      # Converting each dead tree (found in fates_res^2 area) to no of trees per ha that corresponds to the units in the restart file
    # saving plant density to remove for each cohort.rowid ()this matches the sequence of cohorts in the restart file, from where the bft was created in extract.treelist.R
    fire.dead.nplant <- aggregate(remove_nplant ~ fates_pft + fates_dbh + cohort.rowid + nsam, 
			      data = fire.dead.plantlist, FUN = sum, na.rm = TRUE)
    write.table(fire.dead.plantlist,
	    file = file.path(outdir, paste0("/fire.dead.plantlist_", cycle_index, ".dat")),
	    row.names = FALSE)

    write.table(fire.dead.nplant,
	    file = file.path(outdir, paste0("/fire.dead.nplant_", cycle_index, ".dat")),
	    row.names = FALSE)

    fire.dead.ls <- split(fire.dead.nplant, fire.dead.nplant$nsam)
    # for each nsam in names(fire.dead.ls), go and remove n.plant in nc
    filter.arr <-
      # read.table(file.path(outdir, filterFile), header = F) # on server # commented out by SXM
      read.table(file.path(paste0(outdir, "/", filterFile)), header = F) # ASXM
    sam.vec <-
      c(sam.start:sam.end)[filter.arr$V1]
    nsam <- length(sam.vec)

    # ******************
    # Remove dead nplants from the restart file (reduce densities)
    # By matching cohort.rowid
    # ******************

    # Remove fire-dead trees in units of no. of trees/ha from the restart file.
    # base lists
    pre_post_fire.df.ls <- rep(list(), nsam)
    for (i in 1:length(fire.dead.ls)) {
      sample <- sam.vec[i]
      casename <- paste0(filebase, ".", sample)
      finalhour ="00000"
      filetag <- paste0("elm.r.", finalyear + 1, "-01-01-", finalhour,".nc") # a restart file at the end of final year has timestamp for the beginning of next year
      filename <-
        paste0(runroot,
               "/",
               casename,
               "/run/",
               casename,
               ".",
               filetag)
      # Save filecopy before update
      filetag_for_copy <- paste0("elm.r.", finalyear + 1, "-01-01-", finalhour, "_before_fire.nc")
      filename_for_copy <- paste0(runroot, "/", casename, "/run/", casename, ".", filetag_for_copy)
      file.copy(filename, filename_for_copy)
      print(filename) 
      # Update file
      nc <- nc_open(filename, write = TRUE)
      pre_fire.df <- setNames(data.frame(matrix(ncol = 3, nrow = 1500)), c("fates_pft", "fates_dbh", "fates_nplant"))
      pre_fire.df$fates_pft <- ncvar_get(nc, "fates_pft")
      pre_fire.df$fates_nplant <- ncvar_get(nc, "fates_nplant")
      pre_fire.df$fates_dbh <- ncvar_get(nc, "fates_dbh")
      pre_fire.df$nsam <- i

      post_fire.df <- pre_fire.df
      fire.dead.nplant.i <- fire.dead.ls[[names(fire.dead.ls)[i]]]
      # the above has nplants to remove for each cohort.rowid
      # Reducing tree densities to account for dead trees here
      post_fire.df$fates_nplant[fire.dead.nplant.i$cohort.rowid] <- pre_fire.df$fates_nplant[fire.dead.nplant.i$cohort.rowid] - fire.dead.nplant.i$remove_nplant
      ncdf4::ncvar_put(nc, "fates_nplant", post_fire.df$fates_nplant)
      nc_close(nc)

      # ******************
      # Compare restart file cohorts before and after fire
      # ******************

      colnames(pre_fire.df) <- paste0("pre.", colnames(pre_fire.df))
      colnames(post_fire.df) <- paste0("post.", colnames(post_fire.df))
      pre_post_fire.df <- cbind(pre_fire.df, post_fire.df)
      pre_post_fire.df$diff_nplant <- pre_post_fire.df$pre.fates_nplant - pre_post_fire.df$post.fates_nplant
      pre_post_fire.df.ls[[i]] <- pre_post_fire.df
    }
    pre_post_fire.df <- do.call(dplyr::bind_rows, pre_post_fire.df.ls)
    write.table(pre_post_fire.df,
		  file = file.path(outdir, paste0("/pre_post_fire.df_last_nsam_", cycle_index, ".dat")),
		  row.names = FALSE)
    # Test substitution for the last file:
    nc <- nc_open(filename, write = F)
    ## value in the file
    file.value <- ncvar_get(nc, "fates_nplant")
    ## value that was to be substituted
    desired.value <- post_fire.df$fates_nplant
    #desired.value <- post_fire_n.plant[fire.dead.nplant.i$cohort.rowid]
    ncdf4::nc_close(nc)
    ## make sure the two do match
    match <- all(file.value[fire.dead.nplant.i$cohort.rowid] == desired.value)

    return(match)
    }
  }
