#
# This script load TC track data and create the annual TC affected area  
# Author: Yi-Ying Chen
# First Date: 2019-05-27
#


fun_read_nc <- function(arg1) {
	  #load  ncdf library
	  #
	  library(ncdf4)
  #arg1: filepath fo the nc file from James multilayer output
  print(paste("arg1: for reading file path ;", arg1))
    # open the read in file and copy the variables to the dataframe for analysis
    input_nc <- nc_open(arg1)
    #str(input <- nc)
   #   print(i)

   result <- list()
  for (i in 1:input_nc$ndims ) {
	   #   print(i)
	      # store each  variable with respect of the dim <- name
	      result[[input_nc$dim[[i]]$name]] <- input_nc$dim[[i]]$vals
    }

  for (i in 2:length(input_nc$var) ) {
    # store each variable with respect of the var <- name
    result[[input_nc$var[[i]]$name]] <- ncvar_get(input_nc,input_nc$var[[i]]$name)
  }

  nc_close(input_nc)
  # show result structure
  print(str(result))
  # export the datatable
  return(result)
  }  

lai <- fun_read_nc(arg1="/work/ychen/Satellite/SPOT-VGT/LAI_WEST_PACIFIC/c_gls_LAI_201312310000_WEST_PACIFIC_VGT_V2.0.1.nc")


#print(paste("fun_LAI_dyn works on date:", lai$time[it.tc], sep=" "))


# get dimensions from grid a for working
#mask_$LANDMASK
#nav_lon <- lai$lon
#nav_lat <- lai$lat
nx <- length(lai$lon)
ny <- length(lai$lat)

dis2center <- array(0, dim=c(nx, ny))

radius <-  c(1.0) 
#get land mask 
land.mask <- lai$LAI
land.mask[ is.na(lai$LAI) == FALSE ] <- 1

#create a inital frequency map only for land mask 
land.frq <- land.mask 
land.frq[ is.na(land.mask == 1) ] <- 0 

#pdf("track_mask.pdf",width=8, height=8)
#png("track_mask.png", width = 1024, height = 1024, units = "px")
#layout(matrix(data=seq(1,4),nrow=2, ncol=2, byrow=TRUE))

# load libraries and  TC_track data from the selected year 
library("plyr")
library("fields")

wrk.yr <-seq(1960,1960,1)
run.yr <- 1 

for ( kk in 1:run.yr) { 
    #png(paste("track_mask_",wrk.yr[kk],".png",sep=""), width = 1024, height = 512, units = "px")
    #layout(matrix(data=seq(1,2),nrow=1, ncol=2, byrow=TRUE))
    print(paste("working on year:",wrk.yr[kk],sep=""))

    fname.list <-  list.files(path=paste("/work/ychen/Satellite/SPOT-VGT/TC_track_WP/bwp",wrk.yr[kk],sep=""),
			      pattern = "*.txt|*.dat" )
    # tc.track as a data.frame, update by wrk.yr[kk]  
    tc.track <- data.frame()

    for ( i in 1:length(fname.list) ) {
        track.tmp <- read.csv(file=paste("/work/ychen/Satellite/SPOT-VGT/TC_track_WP/bwp",wrk.yr[kk],"/",fname.list[i],sep="") ,
			      sep=",", header=FALSE)[,1:9] # only import the 1:9 columns  
        #combine the table
        tc.track <- rbind.fill( tc.track, track.tmp)
    }	
      # assign the column names
      colnames(tc.track) <- c("BASIN" ,"CY", "YYYYMMDDHH", "TECHNUM" , "TECH" , "TAU" , "LatNS" , "LonEW" , "VMAX" )

      # seperate text and numbers 
      tc.track$LatNS <- as.numeric( gsub("[^[:digit:]]","",tc.track$LatNS) )/10.
      tc.track$LonEW <- as.numeric( gsub("[^[:digit:]]","",tc.track$LonEW) )/10.

      # subset the dataframe only for the WP basin
        tc.tack <- subset(tc.track, BASIN =="WP")

      # copy tarack information from land.mask
        track.mask <- array(NA, dim=c(nx,ny) )  

      # make a TC mask by the assumption of 50 km  diameter TC track data

	for ( k in 1:length(tc.track$BASIN) ) { 
         #for ( k in 1:1 ) { 
	
         # find TC center location in the map 
           x.center <- tc.track$LonEW[k]
           y.center <- tc.track$LatNS[k]

         # find west/upper and east/lower limits of an retangular domain  
           x1 <- x.center - radius 
           x2 <- x.center + radius
         
           y1 <- y.center - radius
           y2 <- y.center + radius

         # find the TC center lon, lat index with a fixed diameter 
           x.id <- which( (lai$lon > x1) &  (lai$lon < x2) )
           y.id <- which( (lai$lat > y1) &  (lai$lat < y2) )			 

      # find the distance to the center in the square area  
       #for (ii in x.id ) { 
            for (jj in y.id ) { 
	        # update the new track area  
	          dis2center[x.id,jj] <- sqrt( (lai$lon[x.id] - x.center)**2. + (lai$lat[jj] - y.center)**2. )
                  
                   #  if ( dis2center[x.id,jj]  < radius) track.mask[x.id,jj] <- 1                
                   track.mask[x.id[which(dis2center[x.id,jj]<radius)],jj] <- 1
                   # count the accumulative return freq 
                   land.frq[ x.id[which(dis2center[x.id,jj]<radius)],jj ] <-  land.frq[ x.id[which(dis2center[x.id,jj]<radius)],jj ] + 
		                                                              track.mask[x.id[which(dis2center[x.id,jj]<radius)],jj] 
	    }   
       #}
        #lalo <- array(NA, dim=c(length(x.id),length(y.id) ))
        #for (ix in 1:length(x.id)) {
        #    for (iy in 1:length(y.id)) {
        #    lalo
	#    }		    
        #}
        #dis2center[x.id,y.id] <- sqrt( (lai$lon[x.id] - x.center)**2. + (lai$lat[y.id] - y.center)**2. )
        #track.mask[ dis2cemter[x.id,y.id] < radius]  <- 1                
       # assign pixels to the affect area  
       #if (  (length(x.id)>0) &  (length(y.id)>0) ) {  
       # 
       #track.mask [min(x.id):max(x.id), min(y.id):max(y.id)] <- 1   
       #}
       #track.mask[ dis2center < 0.8  ] <- 1
}# end of k for event loop within a year 

track.mask <- track.mask*land.mask
# plot the figure
#image.plot(lai$LAI[,ny:1], main=paste("Typhoon Radius:", radius, "in degree, Year:", wrk.yr[kk], sep=""))  
#image.plot(track.mask[,ny:1], add=T )

# get the information of subset dataset  to Taiwan  
subset.xid <- which( (lai$lon > 119.901785) & (lai$lon < 122.1251))
subset.yid <- which( (lai$lat > 21.83928) & (lai$lat < 25.4018))
track.mask.tw <- track.mask[subset.xid,subset.yid]

#image.plot(track.mask.tw[,length(subset.yid):1], main = paste("Taiwan Region", wrk.yr[kk],sep="") )
# dev.off()
#

# export the track.mask.tw to R data 
# save( file= paste(wrk.yr[kk],"_track_mask.rda",sep=""), track.mask.tw )

# save( file= paste("/work/ychen/scripts/R/Rscripts/SPOT_VGT/west_pacific/",wrk.yr[kk],"_TC_mask.rda",sep=""), track.mask) 


} # end of kk in year loop  

# export the TC affect area return frequency to R data 
# save( file= paste("TC_return_frq_",wrk.yr[1],"to",wrk.yr[run.yr],".rda",sep=""), land.frq )

# dev.off()











