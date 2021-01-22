#Extract the FAO eco zone to the study domain

library(raster)
library(rgdal)
source("levels.R")

ld_go <- T

if (ld_go) {
#p <- shapefile('/lfs/home/ychen/GIS/eco_zone/eco_zone.shp')
fao.eco.shp <- shapefile('/lfs/home/ychen/GIS/eco_zone/eco_zone.shp')
#print(fao.eco.shp)
#assign the prjection and geo reference
fao.eco.shp <- spTransform(fao.eco.shp, CRS('+proj=longlat +datum=WGS84'))
ext <- extent(c(90, 150,0,60))
#ext <- floor(extent(pgeo))

#convert shp to raster 
fao.eco.rr.5km <- raster(ext, res=0.05)
fao.eco.rr.5km <- rasterize(fao.eco.shp, fao.eco.rr.5km, field='GEZ_CODE')

#fine resolution 
fao.eco.rr.1km <- raster(nrow=6722, ncol=6722, ext, res=60/6722)
fao.eco.rr.1km <- resample(x=fao.eco.rr.5km, y=fao.eco.rr.1km, method="ngb")

#convert fao.eco.rr.1km to array 
arr.eco <- as.array(fao.eco.rr.1km)[,,1]
arr.eco <- t(arr.eco)

}#end if ld_go


# load forest cover type
# ESA land cover map from 2015, the domain shold be also matched the LAI map over WP region
load(paste("/lfs/home/ychen/LAI_STUDY_EAsia/LANDCOVER_DATA/2015.esa.landcover.east.asia.rda",sep=""))
# variable name: esa.lc
# group the types to croplands(type 1), forests(type 2), others(type 3)
# see http://http://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-Ph2-PUGv2_2.0.pdf   # Page-30
# set water(210) to NA
esa.lc[ (esa.lc == 210 ) ] <- NA
esa.lc[ (esa.lc >= 10)  & (esa.lc <= 40 ) ] <- 1
esa.lc[ (esa.lc >= 50)  & (esa.lc <= 120) ] <- 2
#esa.lc[ (esa.lc == 160) | (esa.lc == 170) ] <- 2
esa.lc[ (esa.lc > 2) ]  <- 3
# create the LC mask
lc.mask <- esa.lc
# agricuture mask
lc.agr.mask <- lc.mask
lc.agr.mask[lc.mask!=1] <- NA
lc.agr.mask[lc.mask==1] <- 1
# forest mask
lc.for.mask <- lc.mask
lc.for.mask[lc.mask!=2] <- NA
lc.for.mask[lc.mask==2] <- 1
# other mask
lc.oth.mask <- lc.mask
lc.oth.mask[lc.mask!=3] <- NA
lc.oth.mask[lc.mask==3] <- 1


# load TC occerance map already apply for forest only 

run.name <- c("4D_12_100_1999_2018_.rda", "4D_0_100_1999_2018_.rda", "4D_12_0_1999_2018_.rda",
              "3D_10_80_1999_2018_.rda", "3D_0_80_1999_2018_.rda", "3D_10_0_1999_2018_.rda",
              "2D_8_60_1999_2018_.rda", "2D_0_60_1999_2018_.rda", "2D_8_0_1999_2018_.rda")
 

# load runs



tot.aff.runs <- array(0, dim=c(20,9) )


for (irun in 1:9) {
#load the data for iruns
load( file=run.name[irun] )
#decleare array for counting the occurence
count.1999.to.2018 <- array( 0., dim=c(6722,6722))
tmp.arr <- array(0, dim=c(6722,6722))
#sum up all arraies in ann.land.frq
#sum up TC occurence
library("fields")
for ( it in 1:20 ) {
      print(paste("working on year:",it+1999-1) )
      tc.occ.avg[,,it][ tc.occ.avg[,,it]>=1] <- 1.0

      #calculate the total area
      #skip longitude <100E
      tmp.arr <- tc.occ.avg[,,it]*lc.for.mask
      x.id1 = 1  #90E
      x.id2 = 1200 #100E
      tmp.arr [x.id1:x.id2,] <- NA
      #skip logitude <110E & latitude > 30
      #x.id1 =1; x.id2 = 1200*1.5
      #y.id1 =1; y.id2 = 1200*3.0
      #tmp.arr[x.id1:x.id2,y.id1:y.id2]  <- 0
      #count numbers
      #image.plot(tmp.arr[,6722:1])
      #all events
      tmp.arr[tmp.arr > 0.1] <- 1
      tmp.arr[tmp.arr <= 0.1] <- 0
      #plot(tmp.arr)
      #
      tot.aff.runs[it,irun] <- sum(tmp.arr,na.rm=T)

      print(paste("Annual TC affected area:",tot.aff.runs[it,irun],sep=""))
      count.1999.to.2018 <-  count.1999.to.2018 + tmp.arr
      #image.plot(count.1999.to.2018[,6722:1])
}



#count annual occurance of TC (divided by average event numbers from 1999 to 2018, 29 )
tc.ave.occ <- count.1999.to.2018 / 20
#filter out 10 yrs events
#tc.ave.occ[tc.ave.occ<0.1] <- 0


#create tc mask based on the level of TC occurence
tc.mask <- tc.ave.occ 
#set return frequency below ten years  
tc.mask[tc.ave.occ >= 1 ] <- 1
tc.mask[tc.ave.occ <= 0 ] <- NA
image(tc.mask,main="TC mask")
tot.area <- sum( tc.mask, na.rm=T)
print(paste("Total TC affected area:",tot.area,sep="") )

#get TC affect for different eco zone 
tc.eco <- arr.eco*tc.mask 

#go for histgram
aa <- hist(arr.eco, breaks=seq(0,100,1) )
bb <- hist(arr.eco*lc.for.mask, breaks=seq(0,100,1) )
cc <- hist(tc.eco*lc.for.mask, breaks=seq(0,100,1) )
#
tot.area <- 0.
#print unit in Mha
print ( paste( "Run:", run.name[irun], sep="" ) ) 
for( i in 1:length(aa$count) ) {
    if (aa$count[i]/10000 > 0.01) {  
       print(paste("Code:",aa$breaks[i]+1 ,"Count:", format(aa$count[i]/10000, digits=1,nsmall=1),
                   "Forest Count:",format(bb$count[i]/10000,digits=1,nsmall=1),
                   "TC Forest Count:", format(cc$count[i]/10000,digits=1,nsmall=1), sep=","))
       tot.area <- tot.area + (cc$count[i]/10000)
    } #end of if  
} #end of i loop 

print( paste("Sum of total TC affected area:", format(tot.area, digits=1,nsmall=1)))   


} #end of irun






# FAO Global Ecosystem Zone 
#GEZ_TERM		GEZ_CLASS	GEZ_CODE
#Tropical rainforest		Tar	11
#Tropical moist deciduous forestTAwa	12
#Tropical dry forest     	Tawb	13
#Tropical shrubland           	TBsh	14
#Tropical mountain system       TM	16
#Subtropical humid forest	SCf	21
#Subtropical dry forest		SCs	22
#Subtropical steppe		SBSh	23
#Subtropical steppe		SBWh	24
#Tropical desert       		SBWh	24
#Subtropical mountain system	SM	25
#Temperate oceanic forest       TeDo	31
#Temperate continental forest   TeDc	32
#Temperate steppe		TeBSk	33
#Temperate desert              	TeBWk	34
#Temperate mountain system	TeM	35
#Boreal coniferous forest      	Ba	41
#Boreal tundra woodland		Bb	42
#Boreal mountain system		BM	43
#Polar				P	50
#Water				Water	90
#No data        		n.d.	99


