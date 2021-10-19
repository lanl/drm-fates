#-------------------------------------------------
# Find cases that haven't finished with a given (last) time tag
# Modified by Rutuja Chitra-Tarak
# July 25, 2021
#-------------------------------------------------

function_filter <- function(outdir, runroot, filebase, finaltag, start_n, stop_n){
  id.arr <- c(start_n:stop_n)
  nsam <- length(id.arr)
  pb <- txtProgressBar(min = 0, max = nsam, style = 3)
  
  file.rm<-FALSE
  
  filter.arr<-rep(FALSE,nsam)
  filename.arr<-rep(FALSE,nsam)
  
  for(i in 1:nsam){
    setTxtProgressBar(pb, i)
    casename <- paste(filebase, ".", id.arr[i],sep="")
    filename<-paste(runroot,"/",casename,"/run/", casename,".",finaltag,sep="")
    if(file.exists(filename)){
      filter.arr[i] <-TRUE
      filename.arr[i] <- filename
    }
  }
  miss.arr<-id.arr[filter.arr==FALSE]
  write.table(as.data.frame(cbind(filename.arr)),paste0(outdir,"/Filename.txt"),row.names=F,col.names=F)
  write.table(as.data.frame(cbind(filter.arr)),paste0(outdir,"/Filter.txt"),row.names=F,col.names=F)
  write.table(as.data.frame(cbind(miss.arr)),paste0(outdir,"/Missing.txt"),row.names=F,col.names=F)
  return(length(which(filter.arr==TRUE)))
}
