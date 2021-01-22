#
# This script load and analysis SPOT-VGT LAI product
# Author: Yi-Ying Chen
# First Date: 2018-10-06
# Revision: 2019-02-14 

# load the funtion to read the nc file

fun_LAI_dyn <- function(it.tc = 200, wspd.lev=6, wspd.con=5, wspd.ext=15, cut.off=0.1, lai=NA, mask="dynamic") { 

# load library
library("fields")

source("/data1/home/ychen/R/src_function_ncdf4.R")

#it.tc=324
#wspd.lev=4
#wspd.con=1
#wspd.ext=15
#cut.off=0.05
#mask="dynamic" 
#lai=NA

#lai <- fun_read_nc(arg1="/work/ychen/Satellite/SPOT-VGT/LAI_GLOBE_VGT_tw/LAI_19990110to20180620_GLOBE_VGT_tw.nc")

print(paste("fun_LAI_dyn works on date:", lai$time[it.tc], sep=" "))


# get dimensions from grid a for working
#mask <- arr <- a  <- grid <- a$LANDMASK
#nav_lon <- lai$lon
#nav_lat <- lai$lat
nx <- length(lai$lon)
ny <- length(lai$lat)
nt <- length(lai$time)

# set a TC_date

#TC_date=c("2009-07-30")

# find lai observation values only observed before/after 30days of TC_date

# we use 3 dekade files 
# dekade-1 dek_m1
# dekade-2 dek_m2
# dekade-3 dek_m3

# find pixels for dek_m1 
dek_m1 <- array(NA,dim=c(nx,ny)) 
dek_m2 <- array(NA,dim=c(nx,ny)) 
dek_m3 <- array(NA,dim=c(nx,ny)) 
dek_p1 <- array(NA,dim=c(nx,ny)) 
dek_p2 <- array(NA,dim=c(nx,ny)) 
dek_p3 <- array(NA,dim=c(nx,ny)) 
lai_be <- array(NA,dim=c(nx,ny))
lai_af <- array(NA,dim=c(nx,ny))

wspd_mask <- array(NA,dim=c(nx,ny))
wspd_m1 <- array(NA,dim=c(nx,ny)) # ref
wspd_m2 <- array(NA,dim=c(nx,ny)) # windy
# load the WRF wind speed data
# inorder to get exceed wind speed frequency (from Beaufort wind scale 1 to 11)  
# "exc.wspd.frq" 3d matrix(150*180*3)
# and nav.lon; nav.lat 
#it.tc=530 
# check the first 10days of the wrk.year
ld.first.switch <- TRUE
# get working year
wrk.year <- substr(formatC(lai$time[it.tc],width=4, format="d"), start=1,stop=4)
wrk.mon.day <- substr(lai$time[it.tc], start=5,stop=9)
wrk.mon <- as.integer(substr(lai$time[it.tc], start=5,stop=6))
wrk.date <- lai$time[it.tc]
wrk.10d <- (wrk.mon-1)*3 + as.integer(substr(lai$time[it.tc], start=7,stop=7))
# get first 10 days of month output list for working year
if (wrk.mon.day == "0110") {
       ld.first.switch <- TRUE
   }else{
       ld.first.switch <- FALSE
}    
# if ld.first load wind frequency data
#if ( ld.first.switch ) {

load(file=paste(wrk.year,"_wspd.rda",sep=""))
# select the imformation for making wspd_mask 
nav.lon <- output.ls$nav.lon
nav.lat <- output.ls$nav.lat
wsp.frq <- output.ls$exc.frq[,,wspd.lev]
wsp.10d <- output.ls$ws.10day[,,wrk.10d]

load(file=paste(wrk.year,"_rfpd.rda",sep=""))
raf.10d <- output.ls$rf.10day[,,wrk.10d]
# load regrid funtion to resample the 5km grid to 1km grid 
source("/data1/home/ychen/R/fun_regrid_new.R")

#wspd_mask <- t(fun_regrid(nav.lon=nav.lon[,90],nav.lat=nav.lat[75,],array.in=wsp.frq))
#use 99.9% exceed frequency 
wsp.max.10d <- t( fun_regrid(nav.lon=nav.lon[,90],
			      nav.lat=nav.lat[75,],
			      array.in=wsp.10d,type=1))
raf.acc.10d <- t( fun_regrid( nav.lon=nav.lon[,90],
			      nav.lat=nav.lat[75,],
			      array.in=raf.10d,type=1))

# use maximum  10day wind speed to creat wind mask
# low wind 15/ms
ws.test <- 15.0 
wspd_m1[wsp.max.10d <  ws.test] <- 1.0
wspd_m2[wsp.max.10d >= ws.test] <- 1.0 


# use wind probiability during summer  
#wspd.con = wspd.con
#low wind mask
#wspd_m1[wspd_mask ==  0] <- 1.0
#high wind mask
#wspd_m2[wspd_mask >= wspd.con] <- 1.0

#} #ld.first.swith


# load forest map from raw data or R object file 
ld_go <- FALSE
if (ld_go) {
# load Land cover map 
land_cover <- fun_read_nc(arg1=paste("/work/ychen/ycmeet/SPOT_Classification/finalmap/TWD_1997_TM_Taiwan/six_PFT/used_map/",
				     wrk.year,"_500m_6PFT.nc",sep="")) 
# regrid the land cover map to the lai grid

# create a reference data.frame from 2D array: land_cover$lon, land_cover_lat, land_cover$LU_TYPE
land_type_table <- data.frame(lon=c(land_cover$nav_lon), lat=c(land_cover$nav_lat), type=c(land_cover$LU_TYPE))
# land cover class information 
#LU_TYPE=01 for Forest class
#LU_TYPE=02 for Grassland class
#LU_TYPE=04 for Agriculture land class
#LU_TYPE=08 for Built-up land class
#LU_TYPE=13 for Inland water class
#LU_TYPE=26 for Bare soil class
# remove the NA columns in the dataframe 
library(dplyr)
#remove NA
land_type_table <- dplyr::filter(land_type_table,  !is.na(type))
#remove other class only keep FOREST TYPE : type== 1 
land_type_table <- dplyr::filter(land_type_table,  type==1)
print(paste("land_cover_data_information:", dim(land_type_table), sep=""))
print(paste("regriding 500m forest cover map to a 1km forest_mask by using majority approach, set to 50%.", sep=""))

# regriding
spacing <- 0.01/2  # ~1km
for (i in 1:(nx)) {
     for (j in 1:(ny)) {
           #only do the land points 
	   if ( (is.na(lai$LAI[i,j,1])!=TRUE)  )  {  
	     # use which funtion to find id in the dataframe
             tmp_id <- which(land_type_table$lon >= lai$lon[i]-spacing & land_type_table$lon <= lai$lon[i]+spacing &
		             land_type_table$lat >= lai$lat[j]-spacing & land_type_table$lat <= lai$lat[j]+spacing ) 
             # justify the pixel type by counting forest pixel >=4 (100%) 
             if (sum(!is.na(land_type_table$type[tmp_id])) >= 4 ) {
                forest_mask[i,j] <- 1		   
	     }	   
           } 
     }	     
}

# save the forest_mask to R object data 
save(forest_mask, file = "/work/ychen/scripts/R/Rscripts/SPOT_VGT/tw/forest_mask.rda")
} else {
# load the forest mask from R data object from previous calculation
	if (mask == "static") {
        load(file = "/work/ychen/scripts/R/Rscripts/SPOT_VGT/tw/forest_mask.rda")	
           print(paste("use static mask ~ ",sep=""))
	 
	}else{
          if (wrk.year < 2000) yr.id = 1995  
	  if ((wrk.year >= 2000) & (wrk.year <2005)) yr.id = 2000
          if ((wrk.year >= 2005) & (wrk.year <2010)) yr.id = 2005
          if ((wrk.year >= 2010) & (wrk.year <2015)) yr.id = 2010
          if ((wrk.year >= 2015) & (wrk.year< 2020)) yr.id = 2015	  
          print(paste("use mask from year:",yr.id,"~",sep=""))
	  load(file = paste("/work/ychen/scripts/R/Rscripts/SPOT_VGT/tw/",yr.id,"_forest_mask.rda",sep=""))
	}	
}




#





#it.tc = 330
obs_cont = 10
# LAI Change analysis
# loop for spatil domain
for ( ix in 1:nx) {
	for ( iy in 1:ny) {
	
              # get index for time
	      # test use 200
              # dek_m1
		#be_cont=60
		#af_cont=60
	        be_cont = 60  	
		af_cont = 30
		if ( is.na(lai$LENGTH_BEFORE[ix,iy,it.tc])==FALSE &  is.na(lai$NOBS[ix,iy,it.tc])==FALSE ) {
		if ( (lai$LENGTH_BEFORE[ix,iy,it.tc] <= be_cont) & 
	             (lai$LENGTH_AFTER[ix,iy,it.tc] <=  af_cont) &
		     (lai$NOBS[ix,iy,it.tc] >= obs_cont)  )  {
	        	dek_m1[ix,iy] <- 1
		}else {
	        	dek_m1[ix,iy] <- NA
		}	
		}
	      # dek_m2
	        be_cont = 45	
		af_cont = 45
		if ( is.na(lai$LENGTH_BEFORE[ix,iy,it.tc-1])==FALSE &   is.na(lai$NOBS[ix,iy,it.tc-1])==FALSE ) {
		if ( (lai$LENGTH_BEFORE[ix,iy,it.tc-1] <= be_cont) & 
	             (lai$LENGTH_AFTER[ix,iy,it.tc-1] <=  af_cont) &
		     (lai$NOBS[ix,iy,it.tc-1] >= obs_cont) )  {
	        	dek_m2[ix,iy] <- 1
		}else {
	        	dek_m2[ix,iy] <- NA
		}	
		}
              # dek_m3
	        be_cont = 30  	
		af_cont = 60
		if ( is.na(lai$LENGTH_BEFORE[ix,iy,it.tc-2])==FALSE &  is.na(lai$NOBS[ix,iy,it.tc-2])==FALSE ) {
		if ( (lai$LENGTH_BEFORE[ix,iy,it.tc-2] <= be_cont) & 
	             (lai$LENGTH_AFTER[ix,iy,it.tc-2] <=  af_cont) &
		     (lai$NOBS[ix,iy,it.tc-2] >= obs_cont) )  {
	        	dek_m3[ix,iy] <- 1
		}else {
	        	dek_m3[ix,iy] <- NA

		}	
		}
              # calculate the lai before TC date  
 	        	        	
		lai_be[ix, iy] <- mean( c(dek_m1[ix,iy]*lai$LAI[ix,iy,it.tc-0],
					  dek_m2[ix,iy]*lai$LAI[ix,iy,it.tc-1],
					  dek_m3[ix,iy]*lai$LAI[ix,iy,it.tc-2] ), na.rm=T) 


              # dek_p1
	        be_cont = 30  	
		af_cont = 60
		if ( is.na(lai$LENGTH_BEFORE[ix,iy,it.tc])==FALSE &  is.na(lai$NOBS[ix,iy,it.tc])==FALSE) {
		if ( (lai$LENGTH_BEFORE[ix,iy,it.tc] <= be_cont) & 
	             (lai$LENGTH_AFTER[ix,iy,it.tc] <=  af_cont) &
		     (lai$NOBS[ix,iy,it.tc] >= obs_cont) )  {
	        	dek_p1[ix,iy] <- 1
		}else {
	        	dek_p1[ix,iy] <- NA
		}	
		}
	      # dek_p2
	        be_cont = 45  	
	      	af_cont = 45
		if ( is.na(lai$LENGTH_BEFORE[ix,iy,it.tc+1])==FALSE &  is.na(lai$NOBS[ix,iy,it.tc+1])==FALSE ) {
		if ( (lai$LENGTH_BEFORE[ix,iy,it.tc+1] <= be_cont) & 
	             (lai$LENGTH_AFTER[ix,iy,it.tc+1] <=  af_cont) & 
		     (lai$NOBS[ix,iy,it.tc+1] >= obs_cont)  )  {
	        	dek_p2[ix,iy] <- 1
		}else {
	        	dek_p2[ix,iy] <- NA
		}	
		}
              # dek_p3
	        be_cont = 60  	
	        af_cont = 30
		if ( is.na(lai$LENGTH_BEFORE[ix,iy,it.tc+2])==FALSE &  is.na(lai$NOBS[ix,iy,it.tc+2])==FALSE) {
		if ( (lai$LENGTH_BEFORE[ix,iy,it.tc+2] <= be_cont) & 
	             (lai$LENGTH_AFTER[ix,iy,it.tc+2] <=  af_cont) &
		     (lai$NOBS[ix,iy,it.tc+2] >= obs_cont) )  {
	        	dek_p3[ix,iy] <- 1
		}else {
	        	dek_p3[ix,iy] <- NA

		}	
		}
              # calculate the lai before TC date  
 	        	        	
		lai_af[ix, iy] <- mean( c(dek_p1[ix,iy]*lai$LAI[ix,iy,it.tc+0],
					  dek_p2[ix,iy]*lai$LAI[ix,iy,it.tc+1],
					  dek_p3[ix,iy]*lai$LAI[ix,iy,it.tc+2] ), na.rm=T) 

	}
}

#apply forest map to the lai map
lai_be <- lai_be * forest_mask
lai_af <- lai_af * forest_mask

#cut off valuse was set to default as 0.1 or can be adjust by assigned vaule to the function 
#cut.off <- 0.1

del_lai <- (lai_af - lai_be) 
#get the current year 

# create pdf obj for ploting 

pdf(paste("./pdf_plot/",wrk.date,"_plot.pdf",sep=""),width=8, height=12)

layout(matrix(data=seq(1,4,1), nrow=2, ncol=2 ),
              widths=c(1,1,1), heights=c(1,1,1))
par( oma=c(2,2,2,2), mgp=c(2,1,0)) # for the axis title, axis labels and axis line 


#declear the wp.matrix for analysising the impact of LAI change due to wind and rainfall
# wind speed at Beaufort Scale
# 7.9/S4, 10.7/S5, 13.8/S6, 17.1/S7, 20.7/S8, 24.4/S9 
# rainfal at CWB Scale
# 40/S1,  80/S2,    200/S3,  350/S4,  500/S5, 700/S6  
#wsp.scale <- c(  3.3, 5.5,  7.9, 10.7, 13.8, 17.1, 20.7, 24.4 ) 
#raf.scale <- c(   0,   40,   80,  200,  350,  500,  700, 1500 )

wsp.scale <- c(  5,   6,   8,  10,   12,   13,  14,  15,  16, 17) 
raf.scale <- c( 50, 125, 200, 275,  350,  425, 500, 575, 650, 1000)


#create a table summaries the result of condictional LAI change based on wind and rainfall scales 
table.comb <- data.frame()

for ( wsp.id in 1:9 )  {
   for ( raf.id in 1:9 ) {
      # calculate area with meteorogical condition with LAI observation 
      met_area <- array(NA,dim=c(nx,ny))
      cut_area <- array(NA,dim=c(nx,ny))
      # 
      #united wind and rainf for analysis area -> 99
      met_area[which( ((wsp.scale[wsp.id] > wsp.max.10d ) |  
                      (raf.scale[raf.id] < raf.acc.10d & raf.acc.10d <= raf.scale[raf.id+1])) & abs(del_lai) > cut.off  ) ] <- 99  
      #cut off area
      cut_area[which( ((wsp.scale[wsp.id] > wsp.max.10d ) |  
                      (raf.scale[raf.id] < raf.acc.10d & raf.acc.10d <= raf.scale[raf.id+1])) & abs(del_lai) <= cut.off  ) ] <- 99  
      
      #count the avariable,  cut.off pixels 
      tot.ava.pix <- length( which( met_area == 99)  )   
      tot.cut.pix <- length( which( cut_area == 99) )
      frac.ava.cut <- tot.ava.pix/(tot.ava.pix+tot.cut.pix)*100.0

      #joint area to 3
      met_area[which( (wsp.scale[wsp.id] > wsp.max.10d ) &  
                      (raf.scale[raf.id] < raf.acc.10d & raf.acc.10d <= raf.scale[raf.id+1]) & met_area == 99 )   ] <- 3   
      #wind only ---> 1
      met_area[which( (wsp.scale[wsp.id] > wsp.max.10d )  & met_area == 99  ) ] <- 1 
      #rainf only ---> 2
      met_area[which( (raf.scale[raf.id] < raf.acc.10d  & raf.acc.10d <= raf.scale[raf.id+1])  & met_area == 99  ) ] <- 2 

      # for loop with three different meteorological conditions 
      for (met.id in 1:3 ) { 
      #get total meteorogical condition pixels
      tot.met.pix <- length( which( met_area == met.id))         

      # calculate the met_area meets the loss/gain pixels
      joint.loss.pix <- length( which( met_area == met.id & del_lai < -1*cut.off) )   
      joint.gain.pix <- length( which( met_area == met.id & del_lai > cut.off)    )
      
      mean.loss <- mean( del_lai[which(met_area == met.id & del_lai < -1*cut.off)], na.rm=T)
      mean.gain <- mean( del_lai[which(met_area == met.id & del_lai > cut.off)],    na.rm=T)

      # calculate total gain/loss pixels
      #tot.loss.met.pix <- length( which(del_lai < -1*cut.off & met_area == 1) )   
      #tot.gain.met.pix <- length( which(del_lai > cut.off & met_area == 1 )   ) 
     
      # calculate joint condition with LAI loss/gain pixels  
      frac.loss <- joint.loss.pix/tot.met.pix*100.0
      frac.gain <- joint.gain.pix/tot.met.pix*100.0

      # if fraction id less than 10 percent olny for summer than exit the loop
      if ( tot.met.pix < 24 | is.nan(frac.loss)==TRUE ) break   
   
     
      # summary the information to the table
      table.tmp <- data.frame( date.time=lai$time[it.tc], 
			       wsp.scale=wsp.scale[wsp.id],
			       raf.scale=raf.scale[raf.id],
			       frac.loss=frac.loss, 
			       frac.gain=frac.gain,
			       pix.loss=joint.loss.pix,
			       pix.gain=joint.gain.pix,
			       pix.met=tot.met.pix,
			       pix.ava=tot.ava.pix,
			       pix.cut=tot.cut.pix,
                               frac.ava.cut=frac.ava.cut,
			       mean.loss=mean.loss,
			       mean.gain=mean.gain,
			       met.class=met.id)

      table.comb <- rbind( table.comb, table.tmp)

     } # met.condtion


      # got the range and interval for color bar
      #r.range <- c(minValue(r <- dif), maxValue(r <- dif))
      r.range <- c(1,3) 
      #assign the color palette
      cols<- colorRampPalette(c("forestgreen", "yellow","red" ))(100)
      brks<- seq(r.range[1], r.range[2],length.out = 101)

      image.plot(met_area[,ny:1], legend.only = FALSE,
		 main=paste("Wind Scl.:", wsp.scale[wsp.id],"-",wsp.scale[wsp.id+1],
			     "Rain Acc. <", raf.scale[raf.id]), breaks=brks, col=cols )       
      # overlap rainfall/wind/delta LAI
      contour(raf.acc.10d[,ny:1],
 	add=T, zlim=c(0,1000),
	nlevels=4, col="gray",lwd=.5)

      contour(wsp.max.10d[,ny:1], 
	add=T, zlim=c(0,20),
	nlevels=4,col="black",lwd=.5)

      #contour(del_lai[,ny:1], 
      #  add=T, zlim=c(-2.0,-0.1),
      #	nlevels=1,col="red",lwd=1)
      points(x=(which(met_area[,ny:1] >= 0 & del_lai[,ny:1] < -1*cut.off, arr.ind = T)[,1])/nx,
	     y=(which(met_area[,ny:1] >= 0 & del_lai[,ny:1] < -1*cut.off, arr.ind = T)[,2])/ny,
	     type="p",col="gray",cex=0.1, pch=1,lwd=0.05 )



   } #raf.scale  
} #wsp.scale

#close pdf device
dev.off()


# show the table 

print(table.comb)


# lai change table for positive and negtive pixels
# create a table summaries the result of LAI change based on all condictions  
table2.comb <- data.frame()
#
# count positive pixels  accumulative LAI drop
tot.pix.pos <- length( del_lai[which( (del_lai > cut.off) & is.na(del_lai) != TRUE )  ] )
tot.pix.neg <- length( del_lai[which( (del_lai < -1*cut.off) & is.na(del_lai) != TRUE )  ] )

#initiate variables
acc.lai.drop=0
acc.lai.gain=0

# remove cut-off pixels   
del_lai[ abs(del_lai) <= cut.off ] <- NA

for (ix in 1:nx) {
   for (iy in 1:ny) {
       # check data with value 	   
       if ( is.na(del_lai[ix,iy]) != TRUE) {
	    # calculate acculative lai  
	    if( del_lai[ix,iy] < 0) {
	       acc.lai.drop = acc.lai.drop + del_lai[ix,iy]
            } else {
	       acc.lai.gain = acc.lai.gain + del_lai[ix,iy]
            } # end of del_lai
       } # end of if 
   }# end of iy 
} # end of ix

# output the information to the table
table2.comb <- data.frame(obs.date=lai$time[it.tc],
			 tot.pix.pos=tot.pix.pos,
			 acc.lai.gain=acc.lai.gain,
			 tot.pix.neg=tot.pix.neg,
			 acc.lai.drop=acc.lai.drop)
#
ld_go_test_plot = TRUE
#
if ( ld_go_test_plot) {

library(fields)
	
plot.obj <-  (lai_af - lai_be)
# set a cut-off range for abs = 0.10

plot.obj[ abs(plot.obj) <= 0.1 ] <- NA

# set colors 
colorTable <- designer.colors(21, c( "red","lightgrey", "blue")) 
brks<- c(-2.0,seq( -1.0, -0.1,,10), seq( 0.1, 1.0,,10),2.0)

#image.plot(plot.obj[,ny:1],
#	   breaks=brks, col=colorTable)

#add countour 
#contour(wspd_mask[,ny:1],add=T, nlevels = 10) 
#contour(raf.acc.10d[,ny:1],
#	add=T, zlim=c(0,1000),
#	nlevels=4, col="forestgreen",lwd=1)
#contour(wspd.max.10d[,ny:1], 
#	add=T, zlim=c(0,20),
#	nlevels=4,col="black",lwd=2)

lai_mean <- (lai_af + lai_be)/2.

result <- list(lai.obs.date=c(lai$time[it.tc]), 
	       lai.avg = formatC(mean(lai_mean, na.rm=T), format = "f", digits = 2),
	       lai.be  = formatC(mean(lai_be,   na.rm=T), format = "f", digits = 2),
	       lai.af  = formatC(mean(lai_af,   na.rm=T), format = "f", digits = 2),
	       lai.table = table.comb,
	       lai.table2 = table2.comb  )


}else{	


#library(fields)
par(mfcol=c(3,2),
    plt = c(0.05,0.95,0.05,0.95), 
    oma=c(2,2,2,2), 
    mgp=c(2,1,0), 
    mar=c(2,2,2,2))

colorTable <- designer.colors(11, c( "brown","yellow", "darkgreen")) 
brks<- c(seq( 0, 2.5,,6), seq( 3.5, 6.5,,6))

image.plot(lai$lon,lai$lat[ny:1],lai_be[,ny:1]*forest_mask[,ny:1],
	   breaks=brks, col=colorTable)

plot.obj <-  (lai_af - lai_be)*forest_mask
# set a cut-off range for abs = 0.20

plot.obj[ abs(plot.obj) <= 0.1 ] <- NA

# set colors 
colorTable <- designer.colors(21, c( "red","lightgrey", "blue")) 
brks<- c(-2.0,seq( -1.0, -0.1,,10), seq( 0.1, 1.0,,10),2.0)

image.plot(lai$lon,lai$lat[ny:1],plot.obj[,ny:1],
	   breaks=brks, col=colorTable)



# seperate the data by differnt wind condiction gp1(<=BWS5),gp2(>BWS5), pg3(all)  
del.lai.gp <-  c(mean(plot.obj*wspd_m1, na.rm=T), 
		 mean(plot.obj*wspd_m2, na.rm=T), 
		 mean(plot.obj, na.rm=T) ) 

#filter the data by the threshold reported by the product uncertainty 
#plot.obj[abs(plot.obj) > 6.0] <- NA
#plot.obj[abs(plot.obj) < 0.5] <- NA
image.plot(lai$lon,lai$lat[ny:1], plot.obj[,ny:1]*wspd_m1[,ny:1],
	   breaks=brks, col=colorTable,xlab="low_wind" )

text(x=120.5,y=24.6,labels= paste("del.LAI.gp1.low.wind:",
     formatC(del.lai.gp[1],format="f",digits=2),sep=""),cex=1.0)

image.plot(lai$lon,lai$lat[ny:1], plot.obj[,ny:1]*wspd_m2[,ny:1],
	   breaks=brks, col=colorTable,xlab="high_wind" )
text(x=120.5,y=24.4,labels= paste("del.LAI.gp2.high.wind:",
     formatC(del.lai.gp[2],format="f",digits=2),sep=""),cex=1.0)



print(paste("it.tc_date:",lai$time[it.tc],sep=""))
print(paste("gp1:",formatC(del.lai.gp[1],width=8,digits=2),
	    "gp2:",formatC(del.lai.gp[2],width=8,digits=2),
	    "all:",formatC(del.lai.gp[3],width=8,digits=2),
	    sep=""))

text(x=120.5,y=25,labels= paste(lai$time[it.tc]),cex=1.0)
text(x=120.5,y=24.8,labels= paste("del.LAI.all:",
     formatC(del.lai.gp[3],format="f",digits=2),sep=""),cex=1.0)

#show wind mask m1:low m2:high
image.plot(wspd_m1[,ny:1])
image.plot(wspd_m2[,ny:1])
#
#
#
#send dev to pdf file

#dev.copy2pdf(file = paste(lai$time[it.tc],".pdf",sep=""), height=10, width=15)


result <- list(lai.obs.date=c(lai$time[it.tc]), 
	       lai.gp1=del.lai.gp[1],
	       lai.gp2=del.lai.gp[2],
	       lai.gp3=del.lai.gp[3],
	       lai.xy=plot.obj,
	       lai.table=table.comb)
#par()

} #ld_go_test_plot

return(result) 


} # EndofFunction

# dekade+1 dek_p1
# dekade+2 dek_p2
# dekade+3 dek_p3

# load the climate fields such, 10days_max_wind, 10days_accumulaytive_rainf 

