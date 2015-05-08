#' @keywords internal
lag_dist_diff_irreg = function(spdata)
{	
	locs <- spdata[,1:2]
	z <-  spdata[,3]
	z <-  z-mean(z)
	n <-  dim(spdata)[1]
	index <-  sort(rep(1:n,n))
	splags <-  c()
	prod <-  c()
	origin.pts <-  matrix(data = NA, nrow = 0, ncol = 2)
	for(i in 1:n)
	{
		cpt.mat <-  cbind(rep(locs[i,1], n), rep(locs[i,2], n))
		origin.pts <-  rbind(origin.pts, cpt.mat)
		lags.d <-  cbind(locs[,1] - locs[i,1] , locs[,2] - locs[i,2] )
		splags <-  rbind(splags, lags.d)
		prod <-  c(prod, (z[i] - z)^2)
	}
	x.coord <- origin.pts[,1]
	y.coord <- origin.pts[,2]
	x.lag <- splags[,1]
	y.lag <- splags[,2]
	xts <- prod
	rv <- cbind(x.coord, y.coord, index, x.lag, y.lag, xts)
	row.names(rv) <- NULL
	return(rv)
}
#' @keywords internal
epkern = function(u)
{
	rv <- rep(0, length(u))
	good <- which(u < 1 & u > -1)
	rv[good] <- 0.75*(1 - u[good]^2)
	return(rv)	
}
#' @keywords internal
normkern = function(u, tr)
{
	if(tr <= 0)
	{
		stop("invalid truncation parameter")
	}
	rv <- rep(0, length(u))
	good <- which(u < tr & u > -tr)
	norm.const <- integrate(dnorm, -tr, tr)
	rv[good] <-  (1/(sqrt(2*pi))*exp(u[good]^2/(-2)))/norm.const$value
	return(rv)
}
#' @keywords internal
coskern = function(u)
{
	rv <- rep(0, length(u))
	good <- which(u < 1 & u > -1)
	rv[good] <- pi/4 * cos(pi*u[good]/2)
	return(rv)
}

unifkern = function(u)
{
	rv <- rep(0, length(u))
	good <- which(u < 1 & u > -1)
	rv[good] <- 0.5
	return(rv)
}
#' @keywords internal
est_gamma_irreg = function(spdata,lagmat, h = 1, kernel = "norm", truncation = 1)
{
	if(dim(spdata)[2] != 3)
	{
		stop("invalid matrix of spatial data")
	}
	if(dim(lagmat)[2] != 2)
	{
		stop("invalid matrix of spatial lags")
	}
	
	npts <-  dim(spdata)[1]
	nlags <-  dim(lagmat)[1]
	rawdata <- lag_dist_diff_irreg(spdata)
	
	gamma.hat <- apply(lagmat, 1, gamma_hat_t, rawdata, h, kernel, truncation)
	rv <- cbind(lagmat,gamma.hat)
	colnames(rv) <- c("lag.x", "lag.y", "gamma.hat")
	return(rv)
}
#' @keywords internal
gamma_hat_t = function(lag, raw_g_data, h, kernel, truncation)
{
	lag.x <-  lag[1]
	lag.y <-  lag[2]
	
	Dx <-  raw_g_data[,4]
	Dy <-  raw_g_data[,5]	
	
	xarg <- (lag.x - Dx)/h
	yarg <- (lag.y - Dy)/h
	
	if(kernel == "ep")
	{
		top <- epkern(xarg)*epkern(yarg)*raw_g_data[,6]
		bot <- epkern(xarg)*epkern(yarg)
	}
	if(kernel == "norm")
	{
		top <- normkern(xarg, truncation)*normkern(yarg, truncation)*raw_g_data[,6]
		bot <- normkern(xarg, truncation)*normkern(yarg, truncation)
	}
	if(kernel == "cos")
	{
		top <- coskern(xarg)*coskern(yarg)*raw_g_data[,6]
		bot <- coskern(xarg)*coskern(yarg)
	}
	if(kernel == "unif")
	{
		top <- unifkern(xarg)*unifkern(yarg)*raw_g_data[,6]
		bot <- unifkern(xarg)*unifkern(yarg)
	}
	
	ghat <- sum(top)/sum(bot)
	return(ghat/2)	
}
#' @keywords internal
est_gamma3 = function(rawdata, lagmat, h, kernel, truncation)
{	
	gamma.hat <- apply(lagmat, 1, gamma_hat_t, rawdata, h, kernel, truncation)
	rv <- cbind(lagmat,gamma.hat)
	colnames(rv) <- c("lag.x", "lag.y", "gamma.hat")
	return(rv)
}
#' @keywords internal
subfield_coords = function(xlims, ylims, blk.dim, grid.spacing = 1)
{
	if(length(xlims) != 2 | length(ylims) != 2)
	{
		stop("invalid x or y limits")
	}
	
	x.locs <-  seq(xlims[1], xlims[2]-blk.dim[1], by = grid.spacing)
	y.locs <-  seq(ylims[1], ylims[2]-blk.dim[2], by = grid.spacing)
	
	window.coords <- expand.grid(y.locs, x.locs)
	window.coords <- cbind(window.coords[,2], window.coords[,1])
	return(window.coords)
}
#' @keywords internal
get_block_data = function(blk.coord, spdata, blk.dim)
{
	min.x <- blk.coord[1]
	max.x <-  blk.coord[1] + blk.dim[1]
	min.y <- blk.coord[2]
	max.y <- blk.coord[2] + blk.dim[2]
	
	good <-  which(spdata[,1] >= min.x & spdata[,1] < max.x & spdata[,2] >= min.y & spdata[,2] < max.y )
	
	return(spdata[good,])
}
#' @keywords internal
est_subblock_irreg = function(spdata, lagmat,xlims, ylims,blk.width, blk.height, grid.spacing, h, kernel = "unif", truncation = 1)
{
	npts <-  dim(spdata)[1]
	tot.vol <- (xlims[2]-xlims[1])*(ylims[2] - ylims[1])
	blk.vol <- blk.width*blk.height
	prop <- blk.vol/tot.vol
	exp.pts.blk <- prop*npts
	if(exp.pts.blk < 4)
	{ warning("number of expected sampling locations/subblock is less than 4")}
	
	nlags <-  dim(lagmat)[1]
	blk.dim <-  c(blk.width, blk.height)
	
	sb.coords <- subfield_coords(xlims, ylims, blk.dim, grid.spacing)
	n.blks1 <- dim(sb.coords)[1]
	sb.coords.list <- list()
	for(i in 1:n.blks1)
	{
		sb.coords.list[[i]] <- sb.coords[i,]
	}
	
	sb.data <- lapply(sb.coords.list, get_block_data, spdata, blk.dim)
	mat.blks <-  unlist( lapply(sb.data,is.matrix) )
	good <-  which(mat.blks == T)
	sb.data <-  sb.data[good]
	blk.size <- do.call(rbind,lapply(sb.data, dim))
	blk.size <- blk.size[,1]
	bad <- which(blk.size < 4)
	if(length(bad) > 0)
	{
		warning("some subblocks discarded due to inadequate data in the block")
		sb.data <- sb.data[-bad]
	}
	
	n.blks2 <- length(sb.data)
	
	sb.rawdata <- lapply(sb.data, lag_dist_diff_irreg)
	sb.ghats <-  lapply(sb.rawdata, est_gamma3, lagmat, h, kernel, truncation)
	n.lags <- dim(lagmat)[1]
	gamma.hat.mat <- matrix(data = NA, nrow = n.blks2, ncol = n.lags)
	has.na <-  c()
	for(i in 1:n.blks2)
	{
		has.na <- c(has.na, sum(1*is.na(sb.ghats[[i]][,3])) )
		gamma.hat.mat[i,] <- sb.ghats[[i]][,3]
	}
	bad <- which(has.na > 0)
	
	if(length(bad) > 0)
	{
		warning("some subblocks discarded due to under-smoothing (small bandwidth, h)")
		gamma.hat.mat <- gamma.hat.mat[-bad,] 
	}
	n.blks2 <- dim(gamma.hat.mat)[1]
	
	rv <- list("n.good.blks" = n.blks2,"n.bad.blks" = n.blks1-n.blks2, "blk.sizes" = blk.size, "gamma.hats" = gamma.hat.mat )
	
	return(rv)
}