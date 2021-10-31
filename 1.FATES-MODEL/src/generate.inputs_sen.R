#-------------------------------------------------
# Generate surface data files with varying parameters 
# Modified by Rutuja Chitra-Tarak
# July 15, 2021
#-------------------------------------------------

#-------------------------------------------------
# Create the files
generate.surface.files <- function(HPU.tab, PARAM_PATH, surf_basefile) {
  #*******************
  ## Load packages----
  #*******************
  # package ncdf4 and tidyverse are needed and are installed wit hthe conda env
 
  #*******************
  ## Process inputs----
  #******************* 
  n.sam <- nrow(HPU.tab)
  
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

  par.names <- colnames(HPU.tab)
  n.par <- length(par.names)
  
  pb <- txtProgressBar(min = 0, max = n.sam, style = 3)
  for (i in 1:n.sam) {
    setTxtProgressBar(pb, i)
    clonefile2 <- paste0(surf.dir, "/", surf_basefile, i, ".nc")
    surf_para <- ncdf4::nc_open(clonefile2, write = T)
    
    for (j in 1:n.par) {
      var.name <- par.names[j]
      if (var.name != "NA") {
        surf.para.values <- ncdf4::ncvar_get(surf_para, var.name)
        surf.para.values <- HPU.tab[i, par.names[j]]
        ncdf4::ncvar_put(surf_para, var.name, surf.para.values)
      } 
    } #j
    ncdf4::nc_close(surf_para)
  } #i
  ## Check whether a file has been converted as desired:
  clonefile.test <- paste0(surf.dir, "/", surf_basefile, 1, ".nc")
  para.test <- ncdf4::nc_open(clonefile.test, write = T)
  ## value in the file
  file.value <- ncdf4::ncvar_get(para.test, par.names[1])
  ## value that was to be substituted
  desired.value <- HPU.tab[1, par.names[1]]
  ncdf4::nc_close(para.test)
  ## make sure the two do match
  
  match <- all(file.value == desired.value)
  return(match)
}

