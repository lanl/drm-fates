# Generate surface data files with varying parameters 
#-------------------------------------------------

#-------------------------------------------------
# Create the files
generate.surface.files <- function(HPU.tab, path.dir, surf_basefile_rel_path) {
  #*******************
  ## Load packages----
  #*******************
  # package ncdf4 and tidyverse are needed and are installed wit hthe conda env
 
  #*******************
  ## Process inputs----
  #******************* 
  # HPU.tab is expected to have headers: SANDPCT SILTPCT ORGC  
  n.sam <- nrow(HPU.tab)
  
  #*******************
  ## Create clones of surface data file to hold changed parameter----
  #******************* 
  surf.dir <- paste0(path.dir, "/data-raw/surface.sam")
  unlink(surf.dir, recursive=TRUE)
  if(!dir.exists(surf.dir)) {dir.create(surf.dir)}

  for (i in 1:n.sam) {
    filename1 <- paste0(path.dir,"/", surf_basefile_rel_path)
    filename2 <- paste0(surf.dir, "/surfdata_bci_panama_v1_c171113.", i, ".nc")
    file.copy(filename1, filename2, overwrite = T)
  }
  
  #-------------------------------------------------
  #Change the parameter value for texture and organic matter
  #-------------------------------------------------
  
  pb <- txtProgressBar(min = 0, max = n.sam, style = 3)
  for (i in 1:n.sam) {
    setTxtProgressBar(pb, i)
    filename2 <- paste0(surf.dir, "/surfdata_bci_panama_v1_c171113.", i, ".nc")
    surf_para <- ncdf4::nc_open(filename1, write = T)
    ## Create new data to replace those in the default file
    txr.4 <- data.frame(depth = c(0.007, 0.03, 0.062, 0.11, 0.21, 0.36, 0.62, 1.04, 1.73, 2.86),
                        PCT_SAND = rep(HPU.tab$SANDPCT[i], 10),
                        PCT_CLAY = rep(c(100 - c(HPU.tab$SANDPCT[i] + HPU.tab$SILTPCT[i])), 10), 
                        ORGANIC = rep(HPU.tab$ORGC[i], 10))
    txr.par.names <- colnames(txr.4)[-1]
    txr.n.par <- length(txr.par.names)
    
    
    for (j in 1:txr.n.par) {
      var.name <- txr.par.names[j]
      if (var.name != "NA") {
        surf.para.values <- ncdf4::ncvar_get(surf_para, var.name)
        surf.para.values <- txr.4[, txr.par.names[j]]
        ncdf4::ncvar_put(surf_para, var.name, surf.para.values)
      } 
    } #j
    ncdf4::nc_close(surf_para)
  } #i
  ## Check whether a file has been converted as desired:
  filename.test.4 <- paste0(surf.dir, "/surfdata_bci_panama_v1_c171113.", 1, ".nc")
  para.test.4 <- ncdf4::nc_open(filename.test.4, write = T)
  ## value in the file
  file.value.4 <- ncdf4::ncvar_get(para.test.4, txr.par.names[1])
  ## value that was to be substituted
  desired.value.4 <- txr.4[, txr.par.names[1]]
  ncdf4::nc_close(para.test.4)
  ## make sure the two does match
  
  df <- all(file.value.4 == desired.value.4)
  return(df)
}

