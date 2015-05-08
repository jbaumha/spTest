#' Nonparametric Test of Isotropy Using the Sample Semivariogram
#'
#' This function performs the nonparametric test of isotropy from Guan et. al. (2004) for spatial data with uniformly distributed sampling locations. See Guan et. al. (2004) for more details.
#'
#' @export
#' @keywords external
#'
#' @param spdata	An \eqn{n} by \eqn{3} matrix. The first two columns provide \eqn{(x,y)} spatial coordinates. The third column provides data values at the coordinates.
#' @param lagmat A \eqn{k} by \eqn{2} matrix of spatial lags. Each row corresponds to a lag of the form \eqn{(x.lag, y.lag)} for which the semivariogram value will be estimated.
#' @param A	A \eqn{d} by \eqn{k} contrast matrix. The contrasts correspond to contrasts of the estimated semivariogram at the lags given in 'lagmat'.
#' @param df A scalar indicating the row rank of A. This value gives the degrees of freedom for the asymptotic Chi-sq distribution used to compute the p-value.
#' @param subblock.dims	A vector of length two corresponding to the width and height of the moving windows used to estimate the asymptotic variance-covariance matrix. If window width does not evenly divide the width of the sampling region, some data will be ommited during subsampling, i.e., function does not handle partial windows. Same applies to window height and height of sampling region.
#' @param h A scalar giving the bandwidth for the kernel smoother. The same bandwidth is used for lags in both the x and y directions.
#' @param kernel A string taking one of the following values: "norm", "ep", "cos", or "unif", for the normal, Epanechnikov, cosine, or uniform kernel functions. Defaults to normal.
#' @param truncation A scalar providing the truncation value for the normal density if 'kernel' is given as "norm".
#' @param xlims A vector of length two providing the lower and upper x-limits of the sampling region.
#' @param ylims A vector of length two providing the lower and upper y-limits of the sampling region.
#' @param grid.spacing A scalar indicating the spacing of the underlying grid laid on the sampling region to create moving windows. The same spacing is used in the x and y directions.
#' @param subblock.h A scalar giving the bandwidth used for the kernel smoother when estimating the semivariogram on the moving windows (blocks). It is recommended to be less than 1 to maintain nominal test size.
#' @param sig.est.finite Logical. True provides a finite sample correction in estimating Sigma (see Guan et. al. (2004) Section 4.2.2). False provides the empirical variance-covariance matrix of sample semivariogram values computed via the moving windows.
#'
#' @details This function currently only supports square and rectangular sampling regions and does not currently support partial blocks. For example, suppose the sampling region runs from 0 to 20 in the x-direction and from 0 to 30 in the y-direction and an underlying grid of 1 by 1 is laid over the sampling region. Then an ideal value of subblock.dims would be (2,3) since its entries evenly divide the width (20) and height (30), respectively, of the sampling region. Using the vector (3, 4.5) would imply that some data will not be used in the moving windows since these values would create partial blocks in the sampling region. 
#'
#'To preserve the spatial dependence structure, the spatial blocks should have the same shape (i.e. square or rectangle) and orientation as the entire sampling domain.
#'
#' @return \item{gamma.hat}{A matrix of the spatial lags provided and the semivariogram point estimates at those lags used to construct the test statistic.}
#' \item{sigma.hat}{The estimate of asymptotic variance-covariance matrix, Sigma, used to construct test statistic.} 
#' \item{n.subblocks}{The number of moving windows (blocks) used to estimate Sigma.}
#' \item{test.stat}{The calculated test statistic.}
#' \item{pvalue.finite}{The approximate, finite-sample adjusted p-value computed by using the moving windows (see Guan et. al. (2004), Section 3.3 for details).}
#' \item{pvalue.chisq}{The p-value computed using the asymptotic Chi-sq distribution.}
#' \item{h.sb}{The value of h (bandwidth) used to estimate the semivariogram on the moving windows (blocks).}
#'
#' @references Guan, Y., Sherman, M., & Calvin, J. A. (2004). A nonparametric test for spatial isotropy using subsampling. \emph{Journal of the American Statistical Association}, 99(467), 810-821.
#'
#' @seealso \code{\link{MaityTest}} \code{\link{GuanTestGrid}}
#'
#' @examples
#' library(mvtnorm)
#' set.seed(1)
#' #Sample Size
#' N <- 300
#' #Set parameter values for exponential covariance function
#' sigma.sq <- 1
#' tau.sq <- 0.0
#' phi <- 1/4
#' #Generate sampling locations
#' coords <-  cbind(runif(N,0,16), runif(N,0,16))
#' D <-  as.matrix(dist(coords))
#' R <- sigma.sq * exp(-phi*D)
#' R <- R + diag(tau.sq, nrow = N, ncol = N)
#' #Simulate Gaussian spatial data
#' z <- rmvnorm(1,rep(0,N), R, method = "chol")
#' z <- z - mean(z)
#' z <- t(z)
#' mydata <- cbind(coords, z)
#' mylags = rbind(c(1,0), c(0, 1), c(1, 1), c(-1,1))
#' myA = rbind(c(1, -1, 0 , 0), c(0, 0, 1, -1))
#' mysb.dims = c(4,4)
#' myh = 0.7
#' myh.sb = 0.8
#' my.xlims = c(0, 20)
#' my.ylims = c(0, 20)
#' tr <- GuanTestUnif(mydata, mylags, myA, df = 2, mysb.dims, myh, "norm", 1.5,
#'  my.xlims, my.ylims, 1, myh.sb)
#' tr
#'
#' #Simulate data from anisotropic covariance function
#' aniso.angle <- pi/4
#' aniso.ratio <- 2
#' coordsA <- coords.aniso(coords, c(aniso.angle, aniso.ratio))
#' Da <- as.matrix(dist(coordsA))
#' R <- sigma.sq * exp(-phi*Da)
#' R <- R + diag(tau.sq, nrow = N, ncol = N)
#' z <- rmvnorm(1,rep(0,N), R, method = c("chol"))
#' z <-  z-mean(z)
#' z <- t(z)
#' mydata <- cbind(coords, z)
#' #Run the test on the data generated from an anisotropic covariance function
#' tr <- GuanTestUnif(mydata, mylags, myA, df = 2, mysb.dims, myh, "norm", 1.5, 
#' my.xlims, my.ylims, 1, myh.sb)
#' tr

GuanTestUnif = function(spdata, lagmat, A, df, subblock.dims, h = 1, kernel = "norm", truncation = 1.5, xlims, ylims, grid.spacing = 1, subblock.h, sig.est.finite = T)
{
	if(dim(spdata)[2] != 3)
	{stop("matrix of spatial data must have 3 columns")}
	if(dim(spdata)[1] <= 3)
	{stop("matrix of spatial data must have at least 4 rows")}
	if(dim(lagmat)[2] != 2)
	{stop("matrix of spatial lags must have 2 columns")}
	if(dim(lagmat)[1] != dim(A)[2])
	{stop("non-conformable `A` and `lagmat` matrix")}
	if(df <= 0)
	{stop("df must be greater than 0")}
	if(length(subblock.dims) != 2)
	{stop("subblock.dims must be length 2")}
	if(subblock.dims[1] <= 0 | subblock.dims[2] <= 0)
	{stop("subblock dimensions must be positive")}
	if(xlims[1] >= xlims[2])
	{stop("invalid x limits of sampling region")}
	if(ylims[1] >= ylims[2])
	{stop("invalid y limits of sampling region")}
	if(df <= 0)
	{stop("df must be greater than 0")}
	if(grid.spacing <= 0)
	{stop("grid.spacing value must be greater than 0")}
	if(kernel != "norm" & kernel != "ep" & kernel != "cos" & kernel != "unif")
	{stop("invalid kernel name entered")}
	if(truncation <= 0)
	{stop("trunction parameter must be positive")}
	if(h <= 0)
	{stop("bandwidth must be positve")}
	if(subblock.h <= 0)
	{stop("subblock bandwidth must be positive")}
	if(subblock.h > 1)
	{warning("Recommend subblock.h to be less than 1 to maintain nominal size")}
	
	
	if( (xlims[2]-xlims[1]) <= subblock.dims[1])
	{stop("subblock width must be less than the width of sampling region")}
	if( (ylims[2]-ylims[1]) <= subblock.dims[2])
	{stop("subblock height must be less than the height of sampling region")}

	if( ((xlims[2]-xlims[1])%%subblock.dims[1]) != 0 )
	{warning("width of subblock does not divide width of sampling region evenly, some data will be ommited during subsampling (function does not handle incomplete subbocks)")}
	if( ((ylims[2]-ylims[1])%%subblock.dims[2]) != 0 )
	{warning("height of subblock does not divide height of sampling region evenly, some data will be ommited during subsampling (function does not handle incomplete subbocks)")}
	
	n.lags <- dim(lagmat)[1]
	
	rawdata <- lag_dist_diff_irreg(spdata)
	gamma.hat <- est_gamma3(rawdata, lagmat, h, kernel, truncation)
	gh <- gamma.hat[,3]
	
	sig.data <- est_subblock_irreg(spdata, lagmat, xlims, ylims, subblock.dims[1], subblock.dims[2], grid.spacing, subblock.h, kernel, truncation)
	
	block.ghats <- sig.data$gamma.hats
	ghat.mean <-  apply(block.ghats, 2, mean)
	block.ghats.mc <-  block.ghats
	for(i in 1:dim(block.ghats)[1])
	{
		block.ghats.mc[i,] <-  block.ghats[i,] - ghat.mean
	}
	
	#not recommended
	if(sig.est.finite == F)
	{
		sigma.hat <- cov(block.ghats.mc)
	}
	
	if(sig.est.finite == T)
	{
		vol.Dn <- (xlims[2]-xlims[1])*(ylims[2]-ylims[1])
		vol.Dni <- subblock.dims[1]*subblock.dims[2]
		fn <- 1 - (vol.Dni/vol.Dn)
		kn <-  sig.data$n.good.blks
		sigma.hat <- ( vol.Dni/(kn*fn) )*t(block.ghats.mc) %*% block.ghats.mc
	}
	
	vol.Dn <- (xlims[2]-xlims[1])*(ylims[2]-ylims[1])
	
	npts <- dim(spdata)[1]
	
	test.stat <- vol.Dn*t(A%*%gh) %*% solve(A%*%sigma.hat%*%t(A)) %*% (A%*%gh)
	test.stat <- test.stat[1,1]
	if(subblock.h < 1)
	{
		test.stat <- vol.Dn*subblock.h^2*t(A%*%gh) %*% solve(A%*%sigma.hat%*%t(A)) %*% (A%*%gh)
		test.stat <- test.stat[1,1]
	}
		
	blk.sizes <- sig.data$blk.sizes
	n.blks <- sig.data$n.good.blks
	blk.vol <- subblock.dims[1]*subblock.dims[2]
	block.test.stats <-  c()
	for(i in 1:n.blks)
	{
		ghb <- matrix( block.ghats[i,], nrow = n.lags, ncol = 1)
		ts <- blk.sizes[i]*t(A%*%ghb) %*% solve(A%*%sigma.hat%*%t(A)) %*% (A%*%ghb)

		if(subblock.h < 1)
		{
			ts <- blk.vol*subblock.h^2*t(A%*%ghb) %*% solve(A%*%sigma.hat%*%t(A)) %*% (A%*%ghb)
		}
		block.test.stats = c(block.test.stats, ts)
	}

	pvalue.finite <- sum(block.test.stats >= test.stat)/n.blks

	pvalue.chisq <- pchisq(test.stat, df, lower.tail = F)
	
	rv <- list("gamma.hat" = gamma.hat, "sigma.hat" = sigma.hat, "n.subblocks" = n.blks, "test.stat" = test.stat,"pvalue.finite" = pvalue.finite, "pvalue.chisq" = pvalue.chisq, "h.sb" = subblock.h)
}
