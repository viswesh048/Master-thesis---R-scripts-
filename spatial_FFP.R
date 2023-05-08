calc_footprint_FFP_climatology <-
  function (zm, z0, umean, h, ol, sigmav, ustar, wind_dir = 0, 
            domain = c(0, 1000, 0, 1000), dx = NULL, dy = NULL, 
            nx = NULL, ny = NULL, r = NULL, rslayer = NULL,  crop = NULL, smooth_data= NULL,
            pulse = NULL, fig = 1) {
    # Flux footprint climatology based on the simple parameterisation FFP 
    # See Kljun, N., P. Calanca, M.W. Rotach, H.P. Schmid, 2015: 
    # Kljun et al. (2015) FFP was modified to accept a spatial roughness parameter 
    # from LIDAR data to calculate the footprint climatology in this script
    # 
    # This function calculates footprints within a fixed physical domain for a series of
    # time steps, rotates footprints into the corresponding wind direction and aggregates 
    # all footprints to a footprint climatology. The percentage of source area is
    # calculated for the footprint climatology.
    # For determining the optimal extent of the domain (large enough to include footprints)
    # use calc_footprint_FFP.m
    #
    # FFP Input
    #    All vectors need to be of equal length (one value for each time step)
    #    zm       = Measurement height above displacement height (i.e. z-d) [m] 
    #               usually a scalar, but can also be a vector. (EC tower data)
    
    #    z0       = Roughness length [m] value matrix derived from the LIDAR data
    
    #    umean    = Vector of mean wind speed at zm [ms-1]. 
    #    Either z0 or umean is required. If both are given, z0 is selected to calculate the footprint (EC tower data)
    
    #    h        = Vector of boundary layer height [m](EC tower data)
    #    ol       = Vector of Obukhov length [m] (EC tower data)
    #    sigmav   = Vector of standard deviation of lateral velocity fluctuations [ms-1] (EC tower data)
    #    ustar    = Vector of friction velocity [ms-1](EC tower data)
    #    wind_dir = Vector of wind direction in degrees (of 360) for rotation of the footprint (EC tower data)   
    

    #    domain       = Domain size as an array of [xmin xmax ymin ymax] [m]
    #                   Footprint will be calculated for a measurement at [0 0 zm] m
    #                   
    #    dx, dy       = Cell size of domain [m]
    #                   Small dx,dy result in higher spatial resolution and higher computing time
    #                   Default is dx = dy = 2 m 
    #    nx, ny       = Number of grid elements in x and y
    #                   Large nx and ny result in higher spatial resolution and higher computing time
    #                   Default is nx = ny = 1000.
    #    r            = Percentage of source area for which to provide contours. This script could only result for source
    #                   contours up to 70%. 
    #                 
    
    #    rslayer      = Calculate footprint even if zm within roughness sublayer: set rslayer = 1
    #                   Note that this only gives a rough estimate of the footprint as the model is not valid within 
    #                   the roughness sublayer. Default is 0 (i.e. no footprint for within RS).
    #                   z0 is needed for estimation of the RS.
    #    pulse        = Display progress of footprint calculations every pulse-th footprint (e.g., "100")
    #    fig          = Plot an example figure of the resulting footprint (on the screen): set fig = 1. 
    #					Default is 0 (i.e. no figure). 
    #
  
    
    require(EBImage)
    
    #--------------------------------------------------------------------
    # Check input variables
    #--------------------------------------------------------------------
    flag_err <- 0
    ind_return <- 0
    output_list <- checkinput_climat(zm, z0, umean, h, ol, sigmav, ustar, 
                                     wind_dir, domain, dx, dy, nx, ny, r=NULL, rslayer = NULL, smooth_data = 0,
                                     crop = 0, pulse = 0, ind_return, flag_err)
    for (v in 1:length(output_list)) assign(names(output_list)[v], 
                                            output_list[[v]])
    
    #--------------------------------------------------------------------
    # Output variables of FFP and their information
    #--------------------------------------------------------------------
    
    #    FFP          = Structure array with footprint climatology data for measurement at [0 0 zm] m
    FFP <- NULL
    
    #    FFP.x_2d     = x-grid of footprint climatology [m]
    #    FFP.y_2d     = y-grid of footprint climatology [m]
    FFP$x_2d <- NaN
    FFP$y_2d <- NaN
    
    #    FFP.fclim_2d = Normalised footprint function values of footprint climatology [m-2]
    FFP$fclim_2d <- NaN
    
    #    FFP.r        = Percentage of footprint as in input, if provided
    FFP$r <- NULL
    
    #    FFP.fr       = Footprint value at r, if r is provided [m-2]
    #    FFP.xr       = x-array for contour line of r, if r is provided [m]
    #    FFP.yr       = y-array for contour line of r, if r is provided [m]
    #                   For array of percentage values, structure entries can be accessed 
    #                   as FFP(1).r, FFP(1).xr, etc.
    FFP$fr <- NULL
    FFP$xr <- NULL
    FFP$yr <- NULL
    
    #    FFP.n        = Number of footprints calculated and included in footprint climatology
    #    flag_err     = 0 if no error, 1 in case of error, 2 if not all contour plots (r%) within specified domain
    #                   3 if single data points had to be removed (outside validity)
    FFP$n <- NULL
    FFP$flag_err <- flag_err
    if (output_list$ind_return) {
      FFP$flag_err <- 1
    }
    
    #--------------------------------------------------------------------
    # Initialize model variables
    #--------------------------------------------------------------------
    else {
      #--------------------------------------------------------------------
      # Constants from Kljun et al. (2015)
      #--------------------------------------------------------------------
      
      a <- 1.4524
      b <- -1.9914
      c <- 1.4622
      d <- 0.1359
      ac <- 2.17
      bc <- 1.66
      cc <- 20
      
      #limit for neutral scaling
      ol_n <- 5000
      #von Karman
      k <- 0.4
      
      #--------------------------------------------------------------------
      # Determine the grid domain of the footprint
      #--------------------------------------------------------------------
      # Define physical domain in cartesian and polar coordinates
      # Cartesian coordinates
      x <- seq(xmin, xmax, (xmax - xmin)/(nx))
      y <- seq(ymin, ymax, (ymax - ymin)/(ny))
      x_2d <- matrix(rep(x, length(y)), nrow = length(y), ncol = length(x), byrow = T)
      y_2d <- matrix(rep(y, length(x)), nrow = length(y), ncol = length(x))
      
      # Polar coordinates
      # Set theta such that North is pointing upwards and angles increase clockwise
      rho <- sqrt(x_2d^2 + y_2d^2)
      theta <- atan2(x_2d, y_2d)
      
      # initialize raster for footprint climatology
      fclim_2d <- 0 * x_2d
      
      #--------------------------------------------------------------------
      # Start loop on the time series of Eddy covariance data
      #--------------------------------------------------------------------
      for (foot_loop in 1:ts_len) {
        if (pulse != 0) {
          if ((foot_loop%%pulse) == 0) {
            print(paste("Calculating footprint", as.character(foot_loop), 
                        "of", as.character(ts_len)))
          }
        }
        
        if (valid[foot_loop]) {
          
          #--------------------------------------------------------------------
          # Create local (within loop) variables of the necessary parameters (wind direction, obukhov length
          # measurement height, roughness length matrix, friction velocity, PBL height, wind speed 
          # and std dev of lateral velocity fluctuations)
          #--------------------------------------------------------------------
          wind_dirl <- wind_dir[foot_loop]
          oll <- ol[foot_loop]
          zml <- zm[foot_loop]
          if (any(!is.na(z0))) {
            z0l <- z0
          }
          else if (!is.na(umean[foot_loop])) {
            z0l <- NA
            umeanl <- umean[foot_loop]
          }
          hl <- h[foot_loop]
          ustarl <- ustar[foot_loop]
          sigmavl <- sigmav[foot_loop]
          
          #--------------------------------------------------------------------
          # Rotate coordinates into wind direction
          #--------------------------------------------------------------------
          wind_dir_rad <- wind_dirl * pi/180
          thetal <- theta - wind_dir_rad
          
          #--------------------------------------------------------------------
          # Create real scale crosswind integrated footprint and dummy for
          # rotated scaled footprint
          #--------------------------------------------------------------------
          f_ci_dummy <- fstar_ci_dummy <- xstar_ci_dummy <- matrix(0, ncol = ncol(x_2d), 
                                                                   nrow = nrow(x_2d))
          px <- matrix(1, ncol = ncol(x_2d), nrow = nrow(x_2d))  
          
        
          if (!any(is.na(z0l)) & (all(z0l> 0))) {
            
            #--------------------------------------------------------------------
            # Calculating wind shear 
            #-------------------------------------------------------------
            if ((oll <= 0) | (oll >= ol_n)) {
              xx <- (1 - 19 * zml/oll)^0.25
              psi_f <- log((1 + xx^2)/2) + 2 * log((1 + 
                                                      xx)/2) - 2 * atan(xx) + pi/2
            }
            else if ((oll > 0) & (oll < ol_n)) {
              psi_f = -5.3 * zml/oll
            }
            #--------------------------------------------------------------------
            # Calculating upwind distance X* and cross wind integrated footprint Fy*.
            # Replacing scalar roughness length value with matrix roughness length from LIDAR data
            # in Eq (7) & (9) from Kljun et al. (2015)
            #--------------------------------------------------------------------
            
            if (any(log(zml/z0l) > psi_f)) {
              xstar_ci_dummy <- rho * cos(thetal)/zml * 
                (1 - (zml/hl))/(log(zml/(z0l)) - psi_f)
              px <- !is.na(xstar_ci_dummy) & (xstar_ci_dummy > d)
              fstar_ci_dummy[px] <- a * (xstar_ci_dummy[px] - 
                                           d)^b * exp(-c/(xstar_ci_dummy[px] - d))
              f_ci_dummy[px] <- fstar_ci_dummy[px]/zml * 
                (1 - (zml/hl))/(log(zml/t(z0l)) - psi_f)
            }
            else {
              flag_err = 3
              valid[foot_loop] = 0
            }
          }
          
          #--------------------------------------------------------------------
          # If roughness length is missing, friction velocity is used to calculate
          # the footprint climatology based on Eq (6) & (8)
          #--------------------------------------------------------------------
          
          else {
            print(1)
            xstar_ci_dummy <- rho * cos(thetal)/zml * (1 - 
                                                         (zml/hl))/(umeanl/ustarl * k)
            px <- !is.na(xstar_ci_dummy) & (xstar_ci_dummy > d)
            fstar_ci_dummy[px] <- a * (xstar_ci_dummy[px] - 
                                         d)^b * exp(-c/(xstar_ci_dummy[px] - d))
            f_ci_dummy[px] <- fstar_ci_dummy[px]/zml * 
              (1 - (zml/hl))/(umeanl/ustarl * k)
          }
          
          #--------------------------------------------------------------------
          # Calculating cross wind dispersion from Eq (18) of Kljun et al. (2015)
          #--------------------------------------------------------------------
          sigystar_dummy <- 0 * x_2d
          sigystar_dummy[px] <- ac * sqrt(bc * (xstar_ci_dummy[px])^2/(1 + 
                                                                         cc * (xstar_ci_dummy[px])))
          
          if (abs(oll) > ol_n) {
            oll <- -1e+06
          }
          if (oll <= 0) {
            scale_const <- 1e-05 * abs(zml/oll)^(-1) + 0.8
          }
          else {
            scale_const <- 1e-05 * abs(zml/oll)^(-1) + 0.55
          }
          scale_const <- min(c(1, scale_const))
          
          #--------------------------------------------------------------------
          # Calculating std dev of the cross wind distance derived from Eq (13) in Kljun et al. (2015)
          #--------------------------------------------------------------------
          sigy_dummy <- 0 * x_2d
          sigy_dummy[px] <- sigystar_dummy[px]/scale_const * 
            zml * sigmavl/ustarl
          sigy_dummy[sigy_dummy < 0] <- NA
          
          #--------------------------------------------------------------------
          # Calculating the 2-D flux footprint from Eq (10) of Kljun et al. (2015)
          #--------------------------------------------------------------------
          f_2d <- 0 * x_2d
          f_2d[px] <- f_ci_dummy[px]/(sqrt(2 * pi) * sigy_dummy[px]) * 
            exp(-(rho[px] * sin(thetal[px]))^2/(2 * sigy_dummy[px]^2))
          
          #--------------------------------------------------------------------
          # Add to footprint climatology raster
          #--------------------------------------------------------------------
          fclim_2d <- fclim_2d + f_2d
        }
        else {
          print(paste("Skipping footprint calculation:", 
                      as.character(foot_loop)))
        }
      }
      
      if (sum(valid) > 0){
        
        fclim_2d <- fclim_2d/sum(valid)
        
        #--------------------------------------------------------------------
        # Derive footprint ellipsoid incorporating R% of the flux, if requested,
        # starting at peak value. (From Kljun FFP model)
        #--------------------------------------------------------------------
        if (!is.null(r[1])) {
          rs <- r  # assigning user input R % source flux 
        }
        else {
          if (crop == 1) {
            rs <- 0.7
          }
          else {
            rs <- NA
          }
        }
        if (!is.na(rs[1]) & (length(fclim_2d) > 1) & (sum(!is.na(fclim_2d)) > 1)) {
          # Calculate integral of fclim_2d starting at peak value until R% are reached
          FFP$r   <- rs * NA
          FFP$fr  <- rs * NA
          ffp_tmp <- NULL
          
          f_sort  <- sort(c(fclim_2d), decreasing = T) # sorting pdf in decreasing order
          f_sort  <- f_sort[!is.na(f_sort)]
          f_cum   <- cumsum(f_sort) * dx * dy  
          
          for (i in 1:length(rs)) {
            f_diff    <- abs(f_cum - rs[i])  
            ind_r     <- which.min(f_diff) 
            fr        <- f_sort[ind_r] ##finding the minimum value which falls within the R % source
            contour_r <- contourLines(x, y, t(fclim_2d), levels = c(fr)) ## deriving the contour lines 
            new_c <- NULL
            #--------------------------------------------------------------------
            # Retrieving the contour coordinates 
            # starting at peak value
            #--------------------------------------------------------------------
            for (nl in 1:length(contour_r)) { ### removing repetitive values for the contour coordinates 
              c_x <- round(contour_r[[nl]]$x * 10)/10
              c_y <- round(contour_r[[nl]]$y * 10)/10
              if (nl == 1) {
                new_c <- unique(cbind(c_x, c_y))  
                new_c <- rbind(new_c, new_c[1, ]) ## ending the contour at the beginning
              }
              else {
                tmp <- unique(cbind(c_x, c_y))
                tmp <- rbind(tmp, tmp[1, ]) ## ending the contour at the beginning
                new_c <- rbind(new_c, c(NA, NA), tmp)
              }
            }
            
            if (!is.na(r[1])) {
              # No output of R% contour lines if outside domain extent		  
              if ((max(new_c[, 2], na.rm = T) >= max(y_2d, 
                                                     na.rm = T)) | (max(new_c[, 1], na.rm = T) >= 
                                                                    max(x_2d, na.rm = T)) | (min(new_c[, 2], 
                                                                                                 na.rm = T) <= min(y_2d, na.rm = T)) | (min(new_c[, 
                                                                                                                                                  1], na.rm = T) <= min(x_2d, na.rm = T))) {
                flag_err <- 2
                FFP$flag_err <- flag_err
                FFP$r[i] <- NA
                FFP$fr[i] <- NA
                FFP$xr[[i]] <- NA
                FFP$yr[[i]] <- NA
              }
              else {
                FFP$r[i] <- rs[i]
                FFP$fr[i] <- fr
                FFP$xr[[i]] <- c(new_c[, 1])
                FFP$yr[[i]] <- c(new_c[, 2])
              }
            }
          }
        }
        
        #--------------------------------------------------------------------
        # Assign the calculated output variables
        #--------------------------------------------------------------------
        FFP$x_2d <- x_2d
        FFP$y_2d <- y_2d
        FFP$fclim_2d <- t(fclim_2d)
        FFP$n <- sum(valid)
        
        
        #--------------------------------------------------------------------
        # Plotting the Footprints and the probability density function
        #--------------------------------------------------------------------
        if (!is.null(fig)) {
          if (fig == 1) {
            if (!is.na(FFP$fclim_2d[1])) {
              x11()
              image.plot(FFP$x_2d[1,], FFP$y_2d[,1], FFP$fclim_2d, xlab='x [m]',ylab='y [m]')
              if (!is.null(FFP$xr[2])) {
                for (i in 1:length(rs)) lines(FFP$xr[[i]], FFP$yr[[i]], type="l", col="red")
              }
            }
          }
          
        }
        
        
      }
    }
    FFP
}



#--------------------------------------------------------------------
# Function checkinput_climat
#--------------------------------------------------------------------
checkinput_climat <-function (zm, z0, umean, h, ol, sigmav, ustar, wind_dir, domain = c(-1000, 
                                                                      1000, -1000, 1000), dx = NULL, dy = NULL, nx = NULL, ny = NULL, 
            r = NULL, rslayer = NULL, smooth_data = NULL, crop = NULL, 
            pulse = NULL, ind_return=0, flag_err=0) {
    flag_err <- 0
    ind_return <- 0
    ts_len <- length(ustar)
    if (length(ustar) != ts_len) {
      print("Input of vectors only")
    }
    if (length(zm) == 1) {
      zm <- rep(zm, length(ustar))
    }
    else if (length(zm) != ts_len){
      print("zm must be either scalar or of same length as other vectors")
      ind_return <- 1
    }
    if ((length(h) != ts_len) | (length(ol) != ts_len) | (length(sigmav) != 
                                                          ts_len) | (length(wind_dir) != ts_len)) {
      print("Input vectors must be of same length")
      ind_return <- 1
    }
    if (min(zm) <= 0) {
      print("zm must be larger than 0")
      ind_return <- 1
    }
    else if (min(h) < 10) {
      print("h must be larger than 10 m")
      ind_return <- 1
    }
    else if (min(sigmav) < 0) {
      print("sig.v must be larger than 0")
      ind_return <- 1
    }
    else if (min(ustar) < 0) {
      print("ustar must be larger than 0")
      ind_return <- 1
    }
    else if (max(wind_dir) > 360) {
      print("(all) wind direction(s) must be <= 360")
      ind_return <- 1
    }
    else if (min(wind_dir) < 0) {
      print("(all) wind direction(s) must be >= 0")
      ind_return <- 1
    }
    if (is.null(r[1])) {
      r <- seq(10, 70, 10)
    }
    if (!is.na(r[1])) {
      if (max(r) > 1) {
        r <- r/100
      }
      if (max(r) > 0.9) {
        print("R must be ,<= 0.9 or <=90#, larger values were removed")
        r <- r[r <= 0.9]
      }
      r <- sort(r)
    }
    if (is.null(rslayer)) {
      rslayer <- 0
    }
    if (is.null(crop)) {
      crop <- 0
    }
    if (is.null(smooth_data)) {
      smooth_data <- 1
    }
    else if ((smooth_data != 0) && (smooth_data != 1)) {
      print("smooth_data must be 0 or 1")
      ind_return <- 1
    }
    if (is.null(pulse)) {
      if (ts_len <= 20) {
        pulse <- 1
      }
      else {
        pulse <- round(ts_len/100)
      }
    }
    if (is.null(domain[1]) | is.na(domain[1])) {
      xmin <- -1000
      xmax <- 1000
      ymin <- -1000
      ymax <- 1000
    }
    else {
      if (length(domain) != 4) {
        print("domain must be an array of four elements [xmin xmax ymin ymax]")
        ind_return <- 1
      }
      min.extent <- min(domain[4] - domain[3], domain[2] - 
                          domain[1])
      if (min.extent < 1) {
        print("domain extent must be larger than 1 m in both x and y")
        ind_return <- 1
      }
      else {
        xmin <- round(domain[1])
        xmax <- round(domain[2])
        ymin <- round(domain[3])
        ymax <- round(domain[4])
      }
    }
    if (!is.null(dx) | !is.null(dy)) {
      if (is.null(dy)) {
        dy <- dx
      }
      else if (is.null(dx)) {
        dx <- dy
      }
      nx <- round((xmax - xmin)/dx)
      ny <- round((ymax - ymin)/dy)
      if (max(c(nx, ny)) > 2000) {
        print("very small dxy - may slow down calculation and cause problems when plotting")
      }
    }
    else {
      if (is.null(nx) & is.null(ny)) {
        nx <- 1000
        ny <- nx
      }
      else if (is.null(ny)) {
        ny <- nx
      }
      else if (is.null(nx)) {
        nx <- ny
      }
      if (((nx%%1) != 0) | ((ny%%1) != 0)) {
        print("nx (and/or ny) are rounded to next integer value(s)")
        nx <- round(nx)
        ny <- round(ny)
      }
      else if ((min(nx) < 10) | (min(ny) < 10)) {
        print("nx and ny extent must be larger than 10 m in both x and y")
        ind_return <- 1
        nx <- NaN
        ny <- NaN
      }
      else if (max(c(nx, ny)) > 2000) {
        print("very large nx or ny - may slow down calculation and cause problems when plotting")
      }
      dx <- (xmax - xmin)/nx
      dy <- (ymax - ymin)/ny
    }
    if (sum(zm > h) > 0) {
      print("zm must be smaller than h")
      ind_return <- 1
    }
    valid <- rep(1, length(ustar))
    valid[(is.na(ustar)) | (ustar < 0.1) | is.na(ol) | is.na(sigmav) | 
            is.na(wind_dir)] <- 0
    
    list(ind_return = ind_return, flag_err = flag_err, valid = valid, 
         ts_len = ts_len, zm = zm, z0 = z0, wind_dir = wind_dir, 
         xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, nx = nx, 
         ny = ny, dx = dx, dy = dy, r = r, smooth_data = smooth_data, 
         crop = crop, pulse = pulse)
  }

