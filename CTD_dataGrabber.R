ctd_location <- "NBP2113_155081_ctd/data/"

vec <- list.files(ctd_location)[grep(".cnv",list.files(ctd_location))]

dnbp.idx <- grep("dnbp2113",vec)
unbp.idx <- grep("unbp2113",vec)
nbp.idx <- (1:length(vec))[!(1:length(vec) %in% c(dnbp.idx,unbp.idx))]

dnbp.loc <- vec[dnbp.idx]
unbp.loc <- vec[unbp.idx]
nbp.loc <- vec[nbp.idx]


library(oce)

folder_vec <- "NBP2113_155081_ctd/data/"
folder_list <- dnbp.loc

nbp1 <- read.ctd(paste(folder_vec,nbp.loc[3],sep = ""))

summary(nbp1)

plot(nbp1)

plot(nbp1@data$fluorescence,-nbp1@data$depth)

for(i in 1:length(folder_list)){
  if(i == 1){
    row_names <- c("MaxDepth","DCM1","DCM_depth","Lat","Long","SurfaceTemp","Salinity","par")
    col_names <- rep(NA,length(folder_vec))
    
    dfOut <- data.frame(matrix(NA,nrow = length(row_names),ncol = length(folder_list)))
    DateTime <- rep(NA,length(row_names))
  }
  
  file_name <- folder_list[i]
  
  ThisFile <- paste(folder_vec,file_name,sep = "")
  
  CTD_vec <- read.ctd(ThisFile)
  
  col_names[i] <- folder_list[i]
  MaxDepth <- max(CTD_vec@data$depth)
  DCM1 <- max(CTD_vec@data$fluorescence[CTD_vec@data$depth > 0],na.rm = T)
  
  DCM_filter <- CTD_vec@data$fluorescence == DCM1
  
  DCM_depth <- CTD_vec@data$depth[DCM_filter & !is.na(CTD_vec@data$fluorescence)]
  if(length(DCM_depth) > 1){
    DCM_depth <- DCM_depth[2]
  }
  
  lat_vec <- round(CTD_vec@metadata$latitude,3)
  long_vec <- round(CTD_vec@metadata$longitude,3)
  temperature <- max(CTD_vec@data$temperature)
  salinity <- max(CTD_vec@data$salinity[CTD_vec@data$salinity<100])
  par <- max(CTD_vec@data$par)
  DateTime[i] <- as.character(CTD_vec@metadata$date)
  
  dfOut[,i] <- c(MaxDepth,
                 DCM1,DCM_depth,
                 lat_vec,long_vec,
                 temperature,
                 salinity,
                 par)
  
  
  # assign(var_name,CTD_vec,envir = globalenv())
  
  if(i == length(folder_list)){
    names(dfOut) <- col_names
    row.names(dfOut) <- row_names
    
    dfOut <- cbind.data.frame(t(dfOut),"DateTime" = DateTime)
  }
}
