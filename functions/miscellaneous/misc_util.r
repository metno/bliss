
#+ check if two rasters match
rasters_match<-function(r1,r2) {
  return ( projection(r1)==projection(r2) &
           extent(r1)==extent(r2) &
           res(r1)[1]==res(r2)[1] & 
           res(r1)[2]==res(r2)[2])
}

#+
set_NAs_to_NULL<-function(x) {
  if (!is.null(x)) {
    if (length(x)>1) {
      ix<-which(is.na(x))
      x[ix]<-NULL
    } else {
      if (is.na(x)) x<-NULL
    }
  }
  x
}

# + replace elements of a string with date-time elements
replaceDate<-function(string=NULL,
                        date.str=NULL,
                        format="%Y-%m-%d %H:%M:%S") {
#------------------------------------------------------------------------------
  if (is.null(string) | is.null(date.str)) return(NULL)
  Rdate<-as.POSIXlt(str2Rdate(date.str,format=format))
  yyyy<-Rdate$year+1900
  mm<-formatC(Rdate$mon+1,width=2,flag="0")
  dd<-formatC(Rdate$mday,width=2,flag="0")
  hh<-formatC(Rdate$hour,width=2,flag="0")
  out<-gsub("yyyy",yyyy,string)
  out<-gsub("mm",formatC(mm,width=2,flag="0"),out)
  out<-gsub("dd",formatC(dd,width=2,flag="0"),out)
  out<-gsub("hh",formatC(hh,width=2,flag="0"),out)
  out
}

#+ find the n-th largest element from each matrix row 
findRow <- function (x,n) {   
# order by row number then by value
  y<-t(x)
  array(y[order(col(y), y)], dim(y))[nrow(y) - (n-1), ]
}

#+ mean radial distance between an observation and its k-th closest obs
dobs_fun<-function(obs,k) {
# NOTE: k=1 will return 0 everywhere (1st closest obs is the obs itself)
  nobs<-length(obs$x)
  if (nobs<k) return(NA)
  # distance matrices
  disth<-(outer(obs$x,obs$x,FUN="-")**2.+
          outer(obs$y,obs$y,FUN="-")**2.)**0.5
  dobsRow<-findRow(x=disth,n=(nobs-k+1))
  mean(dobsRow)
}

#+ mean radial distance between an observation and its k-th closest obs
distance_closest_k_obs<-function(i,
                                 k=10,
                                 sort_decreasing=FALSE) {
#------------------------------------------------------------------------------
# mode == "analysis"
# global variables
#  obs_x/y. observation coordinates
#  grid_x/y. grid point coordinates
#------------------------------------------------------------------------------
  ifelse(length(obs_x)<k,NA,
         sqrt(sort( (obs_x-grid_x[i])**2 + (obs_y-grid_y[i])**2, 
              decreasing = sort_decreasing)[k]) )
}


