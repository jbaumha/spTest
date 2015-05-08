#' Nonparametric Test of Isotropy Using the Sample Covariogram
#'
#' This function performs the nonparametric test of isotropy from Maity and Sherman (2012) for spatial data with sampling locations following any general spatial sampling design. It uses the Epanechnikov kernel function with an empirical bandwidth parameter to smooth over spatial lags. The asymptotic variance-covariance matrix is estimated using the grid based block bootstrap (GBBB) from Lahiri and Zhu (2006). See Maity and Sherman (2012) for more details.
#'
#' @export
#' @keywords external
#'
#' @param spdata	An \eqn{n} by \eqn{3} matrix. The first two columns provide \eqn{(x,y)} spatial coordinates. The third column provides data values at the coordinates.
#' @param lagmat A \eqn{k} by \eqn{2} matrix of spatial lags. Each row corresponds to a lag of the form \eqn{(x.lag, y.lag)} for which the covariogram value will be estimated.
#' @param A	A \eqn{d} by \eqn{k} contrast matrix. The contrasts correspond to contrasts of the estimated covariogram at the lags given in 'lagmat'.
#' @param df A scalar indicating the row rank of A. This value gives the degrees of freedom for the asymptotic Chi-sq distribution used to compute the p-value.
#' @param subblock.dims	A vector of length two corresponding to the width and height of the blocks used in the GBBB. If block width does not evenly divide the width of the sampling region, some data will be ommited during subsampling, i.e., function does not handle partial blocks. Same applies to block height and height of sampling region.
#' @param xlims A vector of length two providing the lower and upper x-limits of the sampling region.
#' @param ylims A vector of length two providing the lower and upper y-limits of the sampling region.
#' @param nBoot A scalar indicating the number of grid based block bootstrap (GBBB) samples to compute (see Lahiri and Zhu (2006) for details on GBBB).
#' @param grid A vector of length two indicating the width and height of the underlying grid laid on the sampling region to create spatial blocks.
#' @param kappa A scalar corresponding to the tuning parameter to adjust empirical bandwidth.
#' @param user.bandwidth Logical. Set to true to manually override the empirical bandwidth.
#' @param bandwidth A vector of length two providing the user-definted bandwidths for smoothing over spatial lags in the x and y directions. Only used when user.bandwidth = TRUE.
#'
#' @details This function currently only supports square and rectangular sampling regions and does not currently support partial blocks. For example, suppose the sampling region runs from 0 to 20 in the x-direction and from 0 to 30 in the y-direction and an underlying grid of 1 by 1 is laid over the sampling region. Then an ideal value of subblock.dims would be (2,3) since its entries evenly divide the width (20) and height (30), respectively, of the sampling region. Using the vector (3, 4.5) would imply that some data will not be used in the GBBB since these values would create partial blocks in the sampling region. To preserve the spatial dependence structure, the spatial blocks should have the same shape (i.e. square or rectangle) and orientation as the entire sampling domain.
#'
#' @return \item{C.hat}{A matrix of the spatial lags provided and the covariogram point estimates at those lags used to construct the test statistic.}
#' \item{V.hat}{The estimate of the asymptotic variance-covariance matrix, V, used to construct the test statistic.}
#' \item{n.boot}{The number of bootstrap samples used to estimate V.}
#' \item{test.stat}{The calculated test statistic.}
#' \item{pvalue.chisq}{The p-value computed using the asymptotic Chi-sq distribution.}
#'
#' @references Maity, A., & Sherman, M. (2012). Testing for spatial isotropy under general designs. \emph{Journal of Statistical Planning and Inference}, 142(5), 1081-1091.
#' @references Lahiri, S. N., & Zhu, J. (2006). Resampling methods for spatial regression models under a class of stochastic designs. \emph{The Annals of Statistics}, 34(4), 1774-1813.
#'
#' @seealso \code{\link{GuanTestUnif}}
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
#' mylags <- rbind(c(1,0), c(0, 1), c(1, 1), c(-1,1))
#' myA <- rbind(c(1, -1, 0 , 0), c(0, 0, 1, -1))
#' mysb.dims <- c(4,4)
#' my.xlims <- c(0, 20)
#' my.ylims <- c(0, 20)
#' #Number of bootstraps for demonstration only.
#' #In practice, use nBoot > 50
#' my.nBoot <- 3
#' tr <- MaityTest(mydata, mylags, myA, df = 2, mysb.dims, my.xlims, my.ylims, 
#' nBoot = my.nBoot, grid = c(1,1))
#' tr
#'
#' ####NOT RUN####
#' # #Simulate data from anisotropic covariance function
#' # aniso.angle <- pi/4
#' # aniso.ratio <- 2
#' # coordsA <- coords.aniso(coords, c(aniso.angle, aniso.ratio))
#' # Da <- as.matrix(dist(coordsA))
#' # R <- sigma.sq * exp(-phi*Da)
#' # R <- R + diag(tau.sq, nrow = N, ncol = N)
#' # z <- rmvnorm(1,rep(0,N), R, method = c("chol"))
#' # z <-  z-mean(z)
#' # z <- t(z)
#' # mydata <- cbind(coords, z)
#' # #Run the test on the data generated from an anisotropic covariance function
#' # tr <- MaityTest(mydata, mylags, myA, df = 2, mysb.dims, my.xlims, my.ylims, 
#' # nBoot = my.nBoot, grid = c(1,1))
#' # tr
MaityTest = function(spdata, lagmat, A, df, subblock.dims, xlims, ylims, nBoot = 100, grid = c(1,1), kappa = 1, user.bandwidth = F, bandwidth = c(1,1))
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
	if(xlims[1] >= xlims[2])
	{stop("invalid x limits of sampling region")}
	if(ylims[1] >= ylims[2])
	{stop("invalid y limits of sampling region")}
	if(nBoot <= 0)
	{stop("invalid number of bootstraps")}
	if(length(subblock.dims) != 2)
	{stop("subblock.dims must be length 2")}
	if(subblock.dims[1] <= 0 | subblock.dims[2] <= 0)
	{stop("subblock dimensions must be positive")}
	if(kappa <= 0)
	{stop("invalid value of kappa")}
	if(subblock.dims[1] < grid[1])
	{stop("blk.dim[1] must be >= grid[1]")}
	if(subblock.dims[2] < grid[2])
	{stop("blk.dim[2] must be >= grid[2]")}
	if(user.bandwidth == T & length(bandwidth != 2))
	{stop("user provided bandwidth must have length 2")}
	
	if(grid[1] > subblock.dims[1])
	{stop("subblock width must be >= underlying grid width")}
	if(grid[2] > subblock.dims[2])
	{stop("subblock height must be >= underlying grid height")}	
	if( (xlims[2]-xlims[1]) <= subblock.dims[1])
	{stop("subblock width must be less than the width of sampling region")}
	if( (ylims[2]-ylims[1]) <= subblock.dims[2])
	{stop("subblock height must be less than the height of sampling region")}

	if( ((xlims[2]-xlims[1])%%subblock.dims[1]) != 0 )
	{warning("width of subblock does not divide width of sampling region evenly, some data will be ommited during subsampling (function does not handle incomplete subbocks)")}
	if( ((ylims[2]-ylims[1])%%subblock.dims[2]) != 0 )
	{warning("height of subblock does not divide height of sampling region evenly, some data will be ommited during subsampling (function does not handle incomplete subbocks)")}
	if(nBoot < 50)
	{warning("at least 50 bootstrap samples are recommended")}

	rawdata <- lag_dist_prod(spdata)
	bad <- which(rawdata[,1] == 0 & rawdata[,2] == 0)
	rawdata <- rawdata[-bad,]
	chat.mat <- est_chat_MS(rawdata, lagmat, kappa = 1, user.bandwidth, bandwidth)
	chat <- chat.mat[,3]
	
	blk.chats <- est_block_chats(lagmat, spdata, nBoot, subblock.dims, xlims, ylims, grid, kappa, user.bandwidth, bandwidth)
	Vhat <- cov(blk.chats)
	Tn <-  t(A %*% chat) %*% solve(A %*% Vhat %*% t(A)) %*% (A %*% chat)
	Tn <- c(Tn)

	pval.chisq <-  pchisq(Tn, df, lower.tail = F)
	pval.chisq <- c(pval.chisq)

	rv <- list("C.hat" = chat.mat, "V.hat" = Vhat, "n.boot" = nBoot, "test.stat" = Tn, "pvalue.chisq" = pval.chisq )
	return(rv)
}



