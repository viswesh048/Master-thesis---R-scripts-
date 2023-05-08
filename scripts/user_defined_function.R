# USER DEFINED FUNCTIONS THAT WERE USED IN THE ANALYSIS

## This function is used to retrieve the separete x and y local coordinates from the footprint climatology and
## combine them
coords_cleaner <- function(model_result){
  longitude_values <- lapply(model_result$xr, lapply, function(x) st_coordinates(ec_tower)[1]+x)
  longitude_values <- lapply(longitude_values, function(x) c(unlist(x)))
  
  latitude_values <- lapply(model_result$yr, lapply, function(x) st_coordinates(ec_tower)[2]+ x)
  latitude_values <- lapply(latitude_values, function(x) c(unlist(x)))
  
  polygon_coords <- list()
  
  for(i in 1:length(latitude_values)){
    polygon_coords[[i]] <- data.frame(longitude_values[i], latitude_values[i])
    colnames(polygon_coords[[i]]) <- c("long", "lat")
  }
  return(polygon_coords)
}

## This function is used to convert the local coordinates to a global coordinates and create polygons of the contours.

polygon_generator <- function(x, buffer){
  #x is a dataframe with coords values and if there is a NA in it that means the table needs to be separated into separate polygons
  if(any(is.na(x$long))){
    split_rows <- which(is.na(x$long))
    if(length(split_rows)==1){
      polygon1 <- x[1:(split_rows-1),]
      polygon2 <- x[(split_rows+1): (length(x$long)),]
      pol1 = st_polygon(list(cbind(polygon1$long, polygon1$lat)) ) %>% st_sfc(crs= st_crs(buffer))
      pol2 = st_polygon(list(cbind(polygon2$long, polygon2$lat)) ) %>% st_sfc(crs= st_crs(buffer))
      if(length(st_contains(pol1, pol2)[[1]])==1){
        complete_pol <- st_difference(pol1, pol2)
      }else{
        complete_pol <- st_union(pol1, pol2)
      }
    }
    else{
      polygon_list_df <- list()
      from_lis <- c(1 , split_rows + 1)
      for(i in 1 : (length(from_lis)-1)){
        polygon_list_df[[i]] <- x[from_lis[i]:(split_rows[i]-1),]
      }
      polygon_list_df[[length(from_lis)]] <- x[from_lis[length(from_lis)]:length(x$long),]
      polygon_list_df <-Filter(function(x) length(x$long) > 3, polygon_list_df)
      polygon_list_df <- lapply(polygon_list_df, function(x) st_polygon(list(cbind(x$long, x$lat)) ) %>% st_sfc(crs= st_crs(buffer)))
      complete_pol <- polygon_list_df[[1]]
      for(i in 2:length(polygon_list_df)){
        if(length(st_contains(complete_pol,  polygon_list_df[[i]])[[1]])==1){
          complete_pol <- st_difference(complete_pol,  polygon_list_df[[i]])
        }else{
          complete_pol <- st_union(complete_pol,  polygon_list_df[[i]])
        }
      }
    }
  }
  else{
    complete_pol <- st_polygon(list(cbind(x$long, x$lat)) ) %>% st_sfc(crs= st_crs(buffer))
  }
  return(complete_pol)
}

## This function creates the datafame of the created polygons
sa_df_generator <- function(sfc_poly){
  sa_polygon_df <- lapply(sfc_poly, st_as_sf)
  for(i in 1:length(sfc_poly)){
    sa_polygon_df[[i]] <- st_as_sf((sfc_poly[[i]]))
    sa_polygon_df[[i]]$sa_perc <- names(sfc_poly[i])
  }
  sa_polygon_df <-lapply(sa_polygon_df, as.data.frame)
  sa_polygon_df <- bind_rows(sa_polygon_df)
  colnames(sa_polygon_df) <- c("geometry", "sa_perc")
  sa_polygon_df <- st_as_sf(sa_polygon_df)
  sa_polygon_df$sa_perc <- as.numeric(sa_polygon_df$sa_perc)
  sa_polygon_df <- sa_polygon_df %>% 
    arrange(-sa_perc)
}


## Converts the probability density function (pdf) matrix to a raster
raster_generator <- function(model_output){
  dat1=list()
  dat1$x <- seq(min(model_output$x_2d), max(model_output$x_2d), 1)
  dat1$y <- seq(min(model_output$y_2d), max(model_output$y_2d), 1)
  dat1$x <- dat1$x + st_coordinates(ec_tower)[1]
  dat1$y <- dat1$y + st_coordinates(ec_tower)[2]
  
  dat1$z <- model_output$fclim_2d
  
  raster_output <- raster(dat1)
  crs(raster_output) <- CRS('+init=EPSG:25832')
  return(raster_output)
}

## This function combines all the above functions to calculate the source area contours and the pdf raster
## for the total, day and night data
total_day_night <- function(total, day, night){
  total_contour <- coords_cleaner(total)
  total_contour <- lapply(total_contour, polygon_generator, buffer)
  names(total_contour) <- seq(10,70, 10)
  total_contour <- sa_df_generator(total_contour)
  total_contour$area <- st_area(total_contour$geometry)
  total_contour$tperiod <- "all"
  
  
  
  day_time_contour <- coords_cleaner(day)
  day_time_contour <- lapply(day_time_contour, polygon_generator, buffer)
  names(day_time_contour) <- seq(10,70, 10)
  day_time_contour <- sa_df_generator(day_time_contour)
  day_time_contour$area <- st_area(day_time_contour$geometry)
  day_time_contour$tperiod <- "day"
  
  
  nig_time_contour <- coords_cleaner(night)
  nig_time_contour <- lapply(nig_time_contour, polygon_generator, buffer)
  names(nig_time_contour) <- seq(10,70, 10)
  nig_time_contour <- sa_df_generator(nig_time_contour)
  nig_time_contour$area <- st_area(nig_time_contour$geometry)
  nig_time_contour$tperiod <- "night"
  
  total_poly <- rbind(total_contour, day_time_contour, nig_time_contour)
  
  total_raster <- raster_generator(total)
  day_raster <- raster_generator(day)
  nig_raster <- raster_generator(night)
  total_list_rast <- list(total_raster, day_raster, nig_raster)
  names(total_list_rast) <- c("all", "day", "night")
  
  return(list(total_poly, total_list_rast))
}



## This function creates the confusion matrix plot
ggplotConfusionMatrix <- function(m){
  mytitle <- paste("Accuracy", percent_format()(m$overall[1]),
                   "Kappa", percent_format()(m$overall[2]))
  p <-
    ggplot(data = as.data.frame(m$table) ,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    geom_text(aes(x = Reference, y = Prediction, label = Freq)) + theme_bw()+
    theme(legend.position = "none") +
    ggtitle(mytitle)+  theme(text = element_text(size = 20), axis.title.x = element_text(size=14, face="bold.italic"),
                             axis.title.y = element_text(size=14, face="bold.italic"))    
  #theme(
  #  plot.title = element_text(size=14, face="bold.italic"),
  #  axis.title.x = element_text(size=14, face="bold"),
  #  axis.title.y = element_text(size=14, face="bold")
  #)
  return(p)
}

## This function is used to extract the microform cover within the footprints

footpr_level <- function(all_masks, grid_image, apc){
  original_masked <- mask(grid_image, all_masks[[1]])
  original_masked <- freq(original_masked, useNA='no')
  original_masked <- as.data.frame(cbind(original_masked, area=round(original_masked[,2] * apc, digits = 1)))
  original_masked$microform <- cover
  original_masked$model <- "original"
  
  dsm_masked <- mask(grid_image, all_masks[[2]])
  dsm_masked <- freq(dsm_masked, useNA='no')
  dsm_masked <- as.data.frame(cbind(dsm_masked, area=round(dsm_masked[,2] * apc, digits = 1)))
  dsm_masked$microform <- cover
  dsm_masked$model <- "dsm"
  
  dtm_masked <- mask(grid_image, all_masks[[3]])
  dtm_masked <- freq(dtm_masked, useNA='no')
  dtm_masked <- as.data.frame(cbind(dtm_masked, area=round(dtm_masked[,2] * apc, digits = 1)))
  dtm_masked$microform <- cover
  dtm_masked$model <- "dtm"
  
  
  total_df <- rbind(original_masked, dsm_masked, dtm_masked)
  return(total_df)
}
## This function is used to estimate the peak value location of the PDF from the tower.
peak_location <- function(raster_img, tower){
  m<- maxValue(raster_img)
  max_value_location <- xyFromCell(raster_img, which(raster_img[]==m))
  max_value_location <- st_point(c(max_value_location[1], max_value_location[2]))
  max_value_location <-st_sfc(max_value_location, crs=st_crs(tower))
  distance <- st_distance(max_value_location, tower)
  return(distance)
  
}