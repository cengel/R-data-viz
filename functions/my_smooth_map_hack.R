# This is a hack of the tmaptools::smooth_map function. 
# I dont want to see the progress bar
# cengel - oct 2017
# from https://rdrr.io/cran/tmaptools/src/R/smooth_map.R
my_smooth_map_hack <- function (shp, var = NULL, nrow = NA, ncol = NA, N = 250000, 
                                unit = "km", unit.size = 1000, smooth.raster = TRUE, nlevels = 5, 
                                style = ifelse(is.null(breaks), "pretty", "fixed"), breaks = NULL, 
                                bandwidth = NA, threshold = 0, cover.type = NA, cover = NULL, 
                                cover.threshold = 0.6, weight = 1, extracting.method = "full", 
                                buffer.width = NA, to.Raster = FALSE) 
{
  is_sf <- inherits(shp, c("sf", "sfc"))
  if (is_sf) 
    shp <- as(shp, "Spatial")
  bbx <- bb(shp)
  prj <- get_projection(shp)
  #pb <- txtProgressBar(style = 3)
  if (!inherits(shp, c("SpatialPoints", "SpatialPolygons", 
                       "SpatialGrid", "Raster"))) {
    stop("shp is not a Raster nor a SpatialPoints, -Polygons, or -Grid object")
  }
  if (!inherits(shp, c("SpatialPoints"))) {
    if (inherits(shp, "Spatial") && !("data" %in% slotNames(shp))) 
      stop("No data found in shape.")
    if (missing(var)) 
      var <- names(shp)[1]
  }
  if (inherits(shp, c("SpatialPoints", "SpatialPolygons"))) {
    bbx <- bb(bbx, ext = -1.05)
    shp@bbox <- bbx
    asp <- get_asp_ratio(shp)
    if (is.na(nrow) || is.na(ncol)) {
      nrow <- round(sqrt(N/asp))
      ncol <- round(N/nrow)
    }
  }
  else {
    if (!inherits(shp, "RasterLayer")) 
      shp <- raster::raster(shp, layer = var)
    ncol <- ncol(shp)
    nrow <- nrow(shp)
  }
  N <- nrow * ncol
  if (!is_projected(shp) || is.na(unit)) {
    warning("shp is not projected; therefore density values cannot be calculated", 
            call. = FALSE)
    cell.area <- 1
  }
  else {
    cell.width <- (bbx[1, 2] - bbx[1, 1])/(unit.size * ncol)
    cell.height <- (bbx[2, 2] - bbx[2, 1])/(unit.size * nrow)
    cell.area <- cell.width * cell.height
  }
  if (is.na(bandwidth[1])) {
    short_side <- min((bbx[, 2] - bbx[, 1])/unit.size)
    bandwidth <- rep(short_side/100, 2)
  }
  else {
    bandwidth <- rep(bandwidth, length.out = 2)
  }
  cover_r <- raster::raster(raster::extent(bbx), nrows = nrow, ncols = ncol, 
                    crs = prj)
  #setTxtProgressBar(pb, 0.1)
  if (is.na(cover.type)) 
    cover.type <- ifelse(inherits(shp, "SpatialPolygons"), 
                         "original", "rect")
  if (missing(cover)) {
    if (cover.type == "rect") {
      cover <- as(raster::extent(bbx), "SpatialPolygons")
      if (!is.na(prj)) 
        cover <- set_projection(cover, current.projection = prj)
      cover_r[] <- TRUE
    }
    else if (cover.type == "original") {
      if (inherits(shp, "Raster")) {
        warning("cover.type=\"original\" only applied to raster output")
        cover <- as(raster::extent(bbx), "SpatialPolygons")
        if (!is.na(prj)) 
          cover <- set_projection(cover, current.projection = prj)
        cover_r <- shp
      }
      else {
        if (inherits(shp, "SpatialPoints")) {
          cover <- gConvexHull(shp)
        }
        else if (inherits(shp, "SpatialPolygons")) {
          cover <- rgeos::gUnaryUnion(shp)
        }
        if (!gIsValid(cover)) 
          cover <- gBuffer(cover, width = 0)
        cover@bbox <- bbx
        cover_r <- poly_to_raster(cover, nrow = nrow, 
                                  ncol = ncol, to.Raster = TRUE)
      }
    }
    else if (cover.type == "smooth") {
      if (!inherits(shp, "Raster")) 
        stop("Raster shape required when cover.type=\"smooth\"")
      cover_list <- smooth_raster_cover(shp, var = var, 
                                        bandwidth = bandwidth * unit.size, threshold = cover.threshold, 
                                        output = c("RasterLayer", "SpatialPolygons"))
      cover_r <- cover_list$RasterLayer
      cover_r[!cover_r[]] <- NA
      cover <- cover_list$SpatialPolygons
    }
  }
  else {
    cover <- rgeos::gUnaryUnion(cover)
    cover_r <- poly_to_raster(cover, nrow = nrow, ncol = ncol, 
                              to.Raster = TRUE)
    bbc <- bb(cover)
    bbx[, 1] <- pmin(bbx[, 1], bbc[, 1])
    bbx[, 2] <- pmin(bbx[, 2], bbc[, 2])
  }
  #setTxtProgressBar(pb, 0.3)
  if (inherits(shp, "SpatialPoints")) {
    co <- sp::coordinates(shp)
    x <- KernSmooth::bkde2D(co, bandwidth = bandwidth * unit.size, gridsize = c(ncol, 
                                                                    nrow), range.x = list(bbx[1, ], bbx[2, ]))
    var <- "count"
  }
  else {
    if (inherits(shp, "SpatialPolygons")) {
      shp@data <- shp@data[, var, drop = FALSE]
      shp <- poly_to_raster(shp, nrow = nrow, ncol = ncol, 
                            copy.data = TRUE, to.Raster = TRUE)
    }
    if (smooth.raster) {
      m <- as.matrix(shp)
      x <- kde2D(m, bandwidth = bandwidth * unit.size, 
                 gridsize = c(ncol, nrow), range.x = list(bbx[1, 
                                                              ], bbx[2, ]))
    }
    else {
      r <- shp
      lvls <- num2breaks(as.vector(r[]), n = nlevels, style = style, 
                         breaks = breaks)$brks
    }
  }
  #setTxtProgressBar(pb, 0.5)
  apply2kde <- inherits(shp, "SpatialPoints") || smooth.raster
  if (apply2kde) {
    r <- raster::raster(raster::extent(bbx), nrows = nrow, ncols = ncol, 
                crs = prj)
    r[] <- as.vector(x$fhat[, ncol(x$fhat):1])
    names(r) <- var
    r[is.na(cover_r[])] <- NA
    if (inherits(shp, "SpatialPoints")) {
      norm_weight <- length(shp) * weight/sum(r[], na.rm = TRUE)
    }
    else {
      norm_weight <- sum(shp[], na.rm = TRUE)/sum(r[], 
                                                  na.rm = TRUE)
    }
    r[] <- r[] * norm_weight/cell.area
    x$fhat <- x$fhat * norm_weight/cell.area
    lvls <- num2breaks(as.vector(x$fhat), n = nlevels, style = style, 
                       breaks = breaks)$brks
    thresLevel <- (lvls[1] == 0 && lvls[2] > threshold && 
                     threshold != 0)
    if (thresLevel) {
      lvls_orig <- lvls
      lvls <- c(lvls[1], threshold, lvls[-1])
    }
    cl <- contourLines(x$x1, x$x2, x$fhat, levels = lvls)
    if (length(cl) < 1L) 
      stop("No iso lines found")
    if (length(cl) > 10000) 
      stop(paste("Number of iso lines over 10000:", length(cl)))
    cl2 <- contour_lines_to_SLDF(cl, proj4string = CRS(prj))
    if (thresLevel) 
      levels(cl2$level) <- c(0, levels(cl2$level)[-1])
  }
  else {
    thresLevel <- FALSE
    bbr <- bb(r)
    rxtra <- (floor(nrow(r)/10) + 1) * 2
    cxtra <- (floor(ncol(r)/10) + 1) * 2
    bbr2 <- bb(bbr, width = (ncol(r) + cxtra)/ncol(r), height = (nrow(r) + 
                                                                   rxtra)/nrow(r), relative = TRUE)
    r2 <- raster::extend(r, raster::extent(bbr2))
    r2[1:(rxtra/2), (cxtra/2 + 1):(ncol(r2) - cxtra/2)] <- r[1, 
                                                             ]
    r2[(nrow(r2) - (rxtra/2) + 1):nrow(r2), (cxtra/2 + 1):(ncol(r2) - 
                                                             cxtra/2)] <- r[nrow(r), ]
    r2[, 1:(cxtra/2)] <- r2[, cxtra/2 + 1]
    r2[, (ncol(r2) - cxtra/2 + 1):ncol(r2)] <- r2[, ncol(r2) - 
                                                    cxtra/2]
    cl2 <- rasterToContour(r2, maxpixels = length(r2), levels = lvls)
  }
  #setTxtProgressBar(pb, 0.7)
  cp <- lines2polygons(ply = cover, lns = cl2, rst = r, lvls = lvls, 
                       extracting.method = "full", buffer.width = buffer.width)
  if (thresLevel) {
    ids <- as.integer(cp$level)
    ids[ids == 1] <- NA
    ids <- ids - 1L
    cp$level <- factor(ids, levels = 1:(length(lvls_orig) - 
                                          1), labels = fancy_breaks(lvls_orig, intervals = TRUE), 
                       ordered = TRUE)
  }
  if (is_sf) 
    cp <- as(cp, "sf")
  attr(cp, "kernel_density") <- TRUE
  #setTxtProgressBar(pb, 0.9)
  lns <- SpatialLinesDataFrame(gIntersection(cover, cl2, byid = TRUE), 
                               data = cl2@data, match.ID = FALSE)
  if (is_sf) 
    lns <- as(lns, "sf")
  attr(lns, "isolines") <- TRUE
  #setTxtProgressBar(pb, 1)
  if (apply2kde && thresLevel) 
    r[][r[] < threshold] <- NA
  list(raster = if (to.Raster) r else as(r, "SpatialGridDataFrame"), 
       iso = lns, polygons = cp, bbox = bbx, nrow = nrow, ncol = ncol, 
       cell.area = cell.area, bandwidth = bandwidth)
}

contour_lines_to_SLDF <- function (cL, proj4string = CRS(as.character(NA)))
{
  require(sp)
  require(rgeos, quietly = TRUE)
  
  .contourLines2LineList <- function (cL)
  {
    n <- length(cL)
    res <- vector(mode = "list", length = n)
    for (i in 1:n) {
      crds <- cbind(cL[[i]][[2]], cL[[i]][[3]])
      res[[i]] <- Line(coords = crds)
    }
    res
  }
  cLstack <- tapply(1:length(cL), sapply(cL, function(x) x[[1]]),
                    function(x) x, simplify = FALSE)
  df <- data.frame(level = factor(names(cLstack), levels=names(cLstack), ordered=TRUE))
  m <- length(cLstack)
  res <- vector(mode = "list", length = m)
  IDs <- paste("C", 1:m, sep = "_")
  row.names(df) <- IDs
  for (i in 1:m) {
    res[[i]] <- Lines(.contourLines2LineList(cL[cLstack[[i]]]),
                      ID = IDs[i])
  }
  SL <- SpatialLines(res, proj4string = proj4string)
  SpatialLinesDataFrame(SL, data = df)
}


buffer_width <- function(bbx) {
  prod(bbx[,2] - bbx[,1]) / 1e12
}


lines2polygons <- function(ply, lns, rst=NULL, lvls, extracting.method="full", buffer.width=NA) {
  prj <- get_projection(ply)
  
  # add a little width to lines
  if (is.na(buffer.width)) buffer.width <- buffer_width(bb(ply))
  suppressWarnings(blpi <- gBuffer(lns, width = buffer.width))
  suppressWarnings(ply <- gBuffer(ply, width = 0))
  
  # cut the poly with isolines
  dpi <- gDifference(ply, blpi)
  
  if (missing(rst)) {
    dpi
  } else {
    # place each polygon in different SpatialPolygon
    ps <- lapply(dpi@polygons[[1]]@Polygons, function(poly) {
      SpatialPolygons(list(Polygons(list(poly), ID = "1")), proj4string = CRS(prj))
    })
    
    # find holes
    holes <- sapply(dpi@polygons[[1]]@Polygons, function(poly) poly@hole)
    
    if (all(holes)) stop("All polygons are holes.")
    
    # create poly id (put each polygon in different feature, and append all holes)
    polyid <- cumsum(!holes)
    
    if (any(holes)) {
      ps_holes <- do.call("sbind", ps[holes])
      ps_solid <- do.call("sbind", ps[!holes])
      
      is_parent <- gContains(ps_solid, ps_holes, byid=TRUE)
      suppressWarnings(areas <- gArea(ps_solid, byid = TRUE))
      parents <- apply(is_parent, MARGIN=1, function(rw) {
        id <- which(rw)
        id[which.min(areas[id])]
      })
      parents <- which(!holes)[parents]
      
      polyid[holes] <- polyid[parents]
    }
    
    m <- max(polyid)
    
    dpi2 <- SpatialPolygons(lapply(1:m, function(i) {
      Polygons(dpi@polygons[[1]]@Polygons[which(polyid==i)], ID=i)
    }), proj4string = CRS(prj))
    
    if (extracting.method=="single") {
      pnts <- gPointOnSurface(dpi2, byid = TRUE)
      values <- raster::extract(rst, pnts)
    } else if (extracting.method=="grid") {
      values <- sapply(1:m, function(i) {
        p <- dpi2[i,]
        rs <- as(raster(extent(p), nrows=10, ncols=10), "SpatialPoints")
        rs@proj4string <- CRS(prj)
        rs <- gIntersection(rs, p)
        if (is.null(rs)) rs <- gPointOnSurface(p) else {
          rs <- sbind(rs, gPointOnSurface(p))
        }
        mean(raster::extract(rst, rs))
      })
    } else {
      # extracting.method=="full"
      values <- sapply(raster::extract(rst, dpi2), function(x)if (is.null(x)) NA else mean(x, na.rm=TRUE))
    }
    
    
    if (length(lvls)==1) {
      lvls <- c(-Inf, lvls, Inf)
    }
    
    # just in case...
    values[is.na(values) | is.nan(values)] <- lvls[1]
    
    brks <- fancy_breaks(lvls, intervals=TRUE)
    
    ids <- cut(values, lvls, include.lowest=TRUE, right=FALSE, labels = FALSE)
    
    if (any(is.na(ids))) stop("raster values not in range")
    if (length(ids)==1) stop("Something went wrong. Probably threshold value too low.")
    
    
    res <- lapply(1:(length(lvls)-1), function(i) {
      if (any(ids==i)) {
        s <- gUnaryUnion(dpi2[ids==i,])
        SpatialPolygonsDataFrame(s, data.frame(level=factor(brks[i], levels=brks, ordered = TRUE)), match.ID = FALSE)
      } else NULL
    })
    res <- res[!sapply(res, is.null)]
    
    x <- do.call("sbind", res)
  }
}

## from https://rdrr.io/cran/tmap/src/R/xxx_scales.R

fancy_breaks <- function(vec, intervals=FALSE, interval.closure="left", fun=NULL, scientific=FALSE, text.separator="to", text.less.than="less than", text.or.more="or more", digits=NA, ...) {
  args <- list(...)
  n <- length(vec)
  
  if (!is.null(fun)) {
    x <- do.call(fun, list(vec))
  } else {
    ### analyse the numeric vector
    vec_fin <- unique(vec[!is.infinite(vec)])
    frm <- gsub(" ", "", sprintf("%20.10f", abs(vec_fin)))
    
    # get width before decimal point
    mag <- max(nchar(frm)-11)
    
    # get number of decimals (which is number of decimals in vec, which is reduced when mag is large)
    ndec <- max(10 - nchar(frm) + nchar(sub("0+$","",frm)))
    if (is.na(digits)) {
      digits <- max(min(ndec, 4-mag), 0)
      
      # add sign to frm
      frm_sign <- paste0(ifelse(vec_fin<0, "-", "+"), frm)
      
      # test if number of digits is sufficient for unique labels
      if (!scientific) {
        while (anyDuplicated(substr(frm_sign, 1, nchar(frm_sign)-10 + digits)) && (digits < 10)) {
          digits <- digits + 1
        }
      }
      
    }
    
    if (!scientific) {
      if (mag>11 || (mag > 9 && all(vec - floor(vec/1e9)*1e9 < 1))) {
        vec <- vec / 1e9
        ext <- " bln"
      } else if (mag > 8 || (mag > 6 && all(vec - floor(vec/1e6)*1e6 < 1))) {
        vec <- vec / 1e6
        ext <- " mln"
      } else {
        ext <- ""
      }
      
      # set default values
      if (!("big.mark" %in% names(args))) args$big.mark <- ","
      if (!("format" %in% names(args))) args$format <- "f"
      if (!("preserve.width" %in% names(args))) args$preserve.width <- "none"
      x <- paste(do.call("formatC", c(list(x=vec, digits=digits), args)), ext, sep="")
    } else {
      if (!("format" %in% names(args))) args$format <- "g"
      x <- do.call("formatC", c(list(x=vec, digits=digits), args))
    }
  }
  
  if (intervals) {
    if (scientific) {
      if (interval.closure=="left") {
        lbls <- paste("[", x[-n], ", ", x[-1], ")", sep="")
        lbls[n-1] <- paste(substr(lbls[n-1], 1, nchar(lbls[n-1])-1), "]", sep="")
      } else {
        lbls <- paste("(", x[-n], ", ", x[-1], "]", sep="")
        lbls[1] <- paste("[", substr(lbls[1], 2, nchar(lbls[1])), sep="")
      }
    } else {
      x[vec==-Inf] <- ""
      lbls <- paste(x[-n], x[-1], sep = paste0(" ", text.separator, " "))
      if (vec[1]==-Inf) lbls[1] <- paste(text.less.than, x[2])
      if (vec[n]==Inf) lbls[n-1] <- paste(x[n-1], text.or.more)
    }
  }
  
  if (intervals) lbls else x
}

num2breaks <- function(x, n, style, breaks, approx=FALSE, interval.closure="left") {
  nobs <- sum(!is.na(x))
  # create intervals and assign colors
  if (style=="fixed") {
    q <- list(var=x,
              brks=breaks)
    attr(q, "style") <- "fixed"
    attr(q, "nobs") <- nobs
    attr(q, "intervalClosure") <- interval.closure
    class(q) <- "classIntervals"
  } else {
    if (nobs==0) stop("Numerical variable only contains missing values.", call.=FALSE)
    if (length(x)==1) stop("Numerical variable only contains one value. Please use a constant value instead, or specify breaks", call. = FALSE)
    q <- suppressWarnings(classIntervals(x, n, style= style, intervalClosure=interval.closure))
  }
  
  if (approx && style != "fixed") {
    if (n >= length(unique(x)) && style=="equal") {
      # to prevent classIntervals to set style to "unique"
      q <- list(var=x, brks=seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=n))
      attr(q, "intervalClosure") <- interval.closure
      class(q) <- "classIntervals"
    } else {
      brks <- q$brks
      
      # to prevent ugly rounded breaks such as -.5, .5, ..., 100.5 for n=101
      qm1 <- suppressWarnings(classIntervals(x, n-1, style= style, intervalClosure=interval.closure))
      brksm1 <- qm1$brks
      qp1 <- suppressWarnings(classIntervals(x, n+1, style= style, intervalClosure=interval.closure))
      brksp1 <- qp1$brks
      if (min(brksm1) > min(brks) && max(brksm1) < max(brks)) {
        q <- qm1
      } else if (min(brksp1) > min(brks) && max(brksp1) < max(brks)) {
        q <- qp1
      }
    }
  }
  q
}
