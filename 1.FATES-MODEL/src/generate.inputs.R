#-------------------------------------------------
# Generate surface data files with varying parameters 
# Modified by Rutuja Chitra-Tarak
# July 15, 2021
#-------------------------------------------------

#-------------------------------------------------
# Create the files
generate.surface.files <- function(param.tab, PARAM_PATH, surf_basefile) {
  #*******************
  ## Load packages----
  #*******************
  # package ncdf4 and tidyverse are needed and are installed wit hthe conda env
 
  #*******************
  ## Process inputs----
  #******************* 
  # param.tab is expected to have headers: SANDPCT SILTPCT ORGC  
  n.sam <- nrow(param.tab)
  
  #*******************
  ## Create clones of surface data file to hold changed parameter----
  #******************* 
  surf.dir <- paste0(PARAM_PATH, "/", "surf.params")
  unlink(surf.dir, recursive=TRUE)
  if(!dir.exists(surf.dir)) {dir.create(surf.dir)}

  for (i in 1:n.sam) {
    basefile <- paste0(PARAM_PATH,"/", surf_basefile, "nc")
    clonefile <- paste0(surf.dir, "/", surf_basefile, i, ".nc")
    file.copy(basefile, clonefile, overwrite = T)
  }
  
  #-------------------------------------------------
  #Change the parameter value for texture and organic matter
  #-------------------------------------------------
  
  pb <- txtProgressBar(min = 0, max = n.sam, style = 3)
  for (i in 1:n.sam) {
    setTxtProgressBar(pb, i)
    clonefile2 <- paste0(surf.dir, "/", surf_basefile, i, ".nc")
    surf_para <- ncdf4::nc_open(clonefile2, write = T)
    ## Create new data to replace those in the default file
    txr.4 <- data.frame(depth = c(0.007, 0.03, 0.062, 0.11, 0.21, 0.36, 0.62, 1.04, 1.73, 2.86),
                        PCT_SAND = rep(param.tab$SANDPCT[i], 10),
                        PCT_CLAY = rep(c(100 - c(param.tab$SANDPCT[i] + param.tab$SILTPCT[i])), 10), 
                        ORGANIC = rep(param.tab$ORGC[i], 10))
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
  clonefile.test <- paste0(surf.dir, "/", surf_basefile, 1, ".nc")
  para.test <- ncdf4::nc_open(clonefile.test, write = T)
  ## value in the file
  file.value <- ncdf4::ncvar_get(para.test, txr.par.names[1])
  ## value that was to be substituted
  desired.value <- txr.4[1, txr.par.names[1]]
  ncdf4::nc_close(para.test)
  ## make sure the two do match

  match <- all(file.value == desired.value)
  return(match)
}

