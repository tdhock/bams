exportData <- function
### Save profiles and annotations as csv.
(data.list=NULL,
### List of profiles and annotations to save to disk as csv.
 out.dir=Sys.getenv("HOME"),
### Directory where to save them.
 quote=FALSE,
### passed to write.table
 row.names=FALSE,
### passed to write.table
 col.names=TRUE,
### passed to write.table
 sep=",",
### passed to write.table
 ...
### additional arguments for write.table
 ){
  if(is.null(data.list)){
    data(neuroblastoma,package="neuroblastoma")
    data.list <- neuroblastoma
  }
  for(N in c("profiles","annotations")){
    fn <- sprintf("%s.csv",N)
    write.table(data.list[[N]],file.path(out.dir,fn),
                quote=quote,
                row.names=row.names,
                col.names=col.names,
                sep=sep,
                ...)
  }
}
