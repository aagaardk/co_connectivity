########################################################################################
# Snow Depth in CO for migration mapping
#---------------------------------------------------------------------------------------
#
# Created by Kevin Aagaard in collaboration with Andy Holland and 
# Eric Bergman
#
# Modified: 
Sys.time()
# By:        Kevin Aagaard
#
#
########################################################################################

# Set directory -----------------------------------------------------------
#Set working directory (should be universal)
setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path)) 
getwd()

# Clear workspace and load packages and libraries -------------------------
rm(list=ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.12")
# BiocManager::install("limma")
x = c("timeDate","lubridate","data.table","plyr","lattice","raster",
      "rgeos","maptools","foreign","xts","limma","sp",
      "PerformanceAnalytics","compiler","fields","rpart","RNCEP",
      "ggmap","rgdal","dplyr","tidyr","tmap","geosphere","tgp",
      "animation","colorspace","colorRamps")
# install.packages(x) #THIS MAY TAKE A WHILE

lapply(x, library, character.only=TRUE)
rm(x)

# Set some global conditions ----------------------------------------------
enableJIT(3)
timeit = proc.time()
options(digits=16)

latlong_proj =
  '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'

aea_proj =
  '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0
+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'

# Read in global data -----------------------------------------------------
NODE_DATA_sp = 
  readOGR(
    ".",
    "co_poly_grid"
  )

NODE_DATA = 
  data.table(
    coordinates(NODE_DATA_sp),
    NODE_DATA_sp@data
  )

setnames(
  NODE_DATA,
  c("V1","V2"),
  c("LongUTM","LatUTM")
)

NODE_DATA[
  ,
  X_INDEX := .GRP,
  by = LongUTM
]

NODE_DATA[
  ,
  Y_INDEX := .GRP,
  by = LatUTM
]

crs(NODE_DATA_sp) = aea_proj

NODE_DATA_sp = 
  sp::spTransform(
    NODE_DATA_sp,
    crs(latlong_proj)
  )

NODE_DATA = cbind(
  NODE_DATA,
  coordinates(NODE_DATA_sp)
)

setnames(
  NODE_DATA,
  c("V1","V2"),
  c("LongDegree","LatDegree")
)

NODE_DATA = 
  data.frame(
    NODE_DATA
  )


#
# Global variables --------------------------------------------------------
num_nodes = nrow(NODE_DATA)
hab_dist = matrix(0,num_nodes,3)

min_year = 2010
max_year = 2019


#
# Download data -----------------------------------------------------------
# # The following code need only be run once to create the files. If the files
# # are lost, this will recreate them. Once they have been saved to the wd(),
# # there is no need to run this.
# 
# SnowDepth.extent =
#   NCEP.gather(
#     variable = 'weasd.sfc',
#     level = 'gaussian',
#     months.minmax =
#       c(
#         1,
#         12
#       ),
#     years.minmax =
#       c(
#         min_year,
#         max_year
#       ),
#     lat.southnorth =
#       c(
#         min(NODE_DATA $LatDegree),
#         max(NODE_DATA $LatDegree)
#       ),
#     lon.westeast =
#       c(
#         min(NODE_DATA $LongDegree),
#         max(NODE_DATA $LongDegree)
#       ),
#     reanalysis2 = FALSE,
#     return.units = TRUE
#   )


#
# Perform conversions -----------------------------------------------------
# Convert snow data from snow-water-equivalence to depth
# Snow density ranges from 10 to 400 in our conditions. We're using 257
# to account for extreme conditions near arctic.
# http://disc.sci.gsfc.nasa.gov/hydrology/data-holdings/parameters/
# snow_water_equivalent.html
density = 257

SnowDepth.extent = SnowDepth.extent / density

# Convert snow depth to centimeters from meters
SnowDepth.extent = SnowDepth.extent * 100


# 
# Resample data onto regular grid -----------------------------------------
# SnowDepth.resamp = array(data=NA,
#                          dim=dim(SnowDepth.extent),
#                          dimnames(SnowDepth.extent))
# 
# pb = winProgressBar(title="Resampling and interpolating onto regular grid",
#                     label="0% done",
#                     min=0,
#                     max=100,
#                     initial=0)
# 
# for(i in 1:dim(SnowDepth.extent)[3]){
#   SnowDepth.reg=interp.loess(x=rep(as.numeric(dimnames(SnowDepth.extent)[[2]]),
#                                    each=length(dimnames(SnowDepth.extent)[[1]])),
#                              y=rep(as.numeric(dimnames(SnowDepth.extent)[[1]]),
#                                    length(dimnames(SnowDepth.extent)[[2]])),
#                              z=as.vector(SnowDepth.extent[,,i]), span=0.75,
#                              gridlen=c(length(dimnames(SnowDepth.extent)[[2]]),
#                                        length(dimnames(SnowDepth.extent)[[1]])))
#   
#   SnowDepth.mat=matrix(data=t(SnowDepth.reg$z),
#                        nrow=length(SnowDepth.reg$y),
#                        ncol=length(SnowDepth.reg$x))
#   SnowDepth.mat=apply(SnowDepth.mat, 2, rev)
#   
#   SnowDepth.resamp[,,i]=SnowDepth.mat
#   
#   info = sprintf("%d%% done", round((i/dim(SnowDepth.extent)[3])*100))
#   setWinProgressBar(pb, i/(dim(SnowDepth.extent)[3])*100, label=info)
# }
# close(pb)


# 
# Convert to data.table via data.frame ------------------------------------
# SnowDepth.df = NCEP.array2df(SnowDepth.resamp)
# SnowDepth.df.cols=colnames(SnowDepth.df)


# 
# Save tables -------------------------------------------------------------
# work_dir = getwd()
# 
# # Snow depth
# filename='/CO_snowdepth_data.txt'
# data_saver=file.path(paste(work_dir,filename,sep=''))
# 
# numrows=10000
# chunksize=floor(nrow(SnowDepth.df)/numrows)
# chunks=0:(chunksize-1)
# chunkseq=chunks*numrows
# pb = winProgressBar(title='Saving Snow Depth Data',
#                     label='0% done',
#                     min=0,
#                     max=100,
#                     initial=0)
# 
# for (i in 1:chunksize){
#   SnowDepth.df_chunk=SnowDepth.df[(1+chunkseq[i]):((1+chunkseq[i])+(numrows-1)),]
#   write.table(SnowDepth.df_chunk,
#               file=data_saver,
#               append=TRUE,
#               row.names=FALSE,
#               col.names=FALSE)
#   
#   info = sprintf('%f%% done', round((i/chunksize)*100))
#   setWinProgressBar(pb, i/(chunksize)*100, label=info)
# }
# close(pb)
# 
# SnowDepth.df_remainder=SnowDepth.df[((chunksize*numrows)+1):nrow(SnowDepth.df),]
# 
# write.table(SnowDepth.df_remainder,
#             file=data_saver,
#             append=TRUE,
#             row.names=FALSE,
#             col.names=FALSE)
# 
#  Convert downloaded data from frames to arrays --------------------------
Snow.dt =
  fread(
    'CO_snowdepth_data.txt'
  )

# Rename columns
setnames(
  Snow.dt,
  names(
    Snow.dt
  ),
  c('date',
    'NCEP_LATITUDE','NCEP_LONGITUDE',
    'daily_snow'
  )
)
Snow.dt[
  ,
  date :=
    substr(
      date,
      1,
      10
    )
]

# Correct negative snow depth values
Snow.dt[
  ,
  daily_snow :=
    ifelse(
      daily_snow < 0,
      0,
      daily_snow
    )
]


# 
# Create averaged data sets -----------------------------------------------
# Snow Depth
Snow.dt[,date:=substr(date,6,10)]
MeanSnow.dt=Snow.dt[,round(mean(daily_snow),digits=3),
                    by=.(date,NCEP_LATITUDE,NCEP_LONGITUDE)]
setnames(MeanSnow.dt,"V1","daily_snow")
MeanSnow.dt=MeanSnow.dt[,.(date,NCEP_LATITUDE,NCEP_LONGITUDE,daily_snow)]

daily_mean_snow=
  dcast(MeanSnow.dt,NCEP_LATITUDE+NCEP_LONGITUDE~date,
        value.var='daily_snow')
setkey(daily_mean_snow,NCEP_LATITUDE,NCEP_LONGITUDE)

# Define nodes based on distance from node centroid to sampling location
NODE_DATA = data.table(NODE_DATA)

NODE_DATA_strpd =
  NODE_DATA[
    ,
    .(
      Y_INDEX,
      X_INDEX,
      LatDegree,
      LongDegree
    )
    ]
setnames(
  NODE_DATA_strpd,
  old =
    colnames(NODE_DATA_strpd),
  new =
    c(
      "Y_INDEX",
      "X_INDEX",
      "ND_LATITUDE",
      "ND_LONGITUDE"
    )
)
NODE_DATA_strpd[
  ,
  ND_LONGITUDE :=
    ND_LONGITUDE + 360
  ]

node_locations =
  NODE_DATA_strpd[
    ,
    .(
      ND_LONGITUDE,
      ND_LATITUDE
    )
    ]
node_locations =
  as.matrix(
    node_locations,
    nrow = nrow(NODE_DATA),
    ncol = 2
  )

# Snow Depth
samp_locations = MeanSnow.dt[,.(NCEP_LONGITUDE,NCEP_LATITUDE)]
samp_locations = unique(samp_locations)
samp_locations=as.matrix(samp_locations,nrow=nrow(MeanSnow.dt),ncol=2)

closest_samp=matrix(NA,nrow=nrow(node_locations),ncol=4)
closest_samp[,1:2]=node_locations
for(i in 1:nrow(node_locations)){
  closest_samp[i,3:4]=
    samp_locations[which.min(distHaversine(node_locations[i,],
                                           samp_locations)),]
}

colnames(closest_samp)=c("ND_LONGITUDE","ND_LATITUDE","NCEP_LONGITUDE",
                         "NCEP_LATITUDE")
closest_samp=data.table(closest_samp)

ND_daily_mean_snow =
  data.table(
    matrix(
      0,
      nrow =
        nrow(
          closest_samp
        ),
      ncol =
        ncol(
          daily_mean_snow
        )
    )
  )

ND_daily_mean_snow = data.table(ND_daily_mean_snow)

ND_daily_mean_snow[,`:=`(V1=closest_samp[,NCEP_LATITUDE],
                          V2=closest_samp[,NCEP_LONGITUDE])]

setnames(ND_daily_mean_snow,names(ND_daily_mean_snow),
         colnames(daily_mean_snow))
setkey(ND_daily_mean_snow,NCEP_LATITUDE,NCEP_LONGITUDE)

ND_daily_mean_snow=
  daily_mean_snow[ND_daily_mean_snow][
    ,1:ncol(daily_mean_snow),with=F]

setkey(closest_samp,ND_LONGITUDE,ND_LATITUDE)
setkey(NODE_DATA_strpd,ND_LONGITUDE,ND_LATITUDE)

full_node_coords=closest_samp[NODE_DATA_strpd]

setkey(ND_daily_mean_snow,NCEP_LONGITUDE,NCEP_LATITUDE)
setkey(full_node_coords,NCEP_LONGITUDE,NCEP_LATITUDE)

ND_daily_mean_snow[,`:=`(NCEP_LATITUDE=NULL,NCEP_LONGITUDE=NULL)]

ND_daily_mean_snow=cbind(full_node_coords[,.(Y_INDEX,X_INDEX)],
                          ND_daily_mean_snow)

ND_daily_mean_snow=as.matrix(ND_daily_mean_snow)


#
# Save files --------------------------------------------------------------
# Snow Depth
write.table(ND_daily_mean_snow,
            paste(getwd(), 
                  "/Colorado_mean_snow_node_info.txt",
                  sep=""), 
            sep="\t", 
            col.names=TRUE, 
            row.names=FALSE)


#