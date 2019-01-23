
#' @title make an example multi-view data set for illustration of MOFA
#' @name makeExampleData
#' @description Function to simulate an example multi-view data set according to the generative model of MOFA.
#' @param n_views number of views
#' @param n_features number of features in each view 
#' @param n_samples number of samples
#' @param n_factors number of factors
#' @param likelihood likelihood for each view, one of "gaussian", "bernoulli", "poisson",
#'  or a character vector of length n_views
#' @return Returns an untrained \code{\link{MOFAmodel}} containing simulated data as training data.
#' @importFrom stats rnorm rbinom rpois
#' @export
#' @examples
#' 
#' # Generate a data set
#' MOFAexample <- makeExampleData()
#' MOFAexample

makeExampleData <- function(n_views=3, n_features=100, n_samples = 50,
                            n_factors = 5, likelihood = "gaussian") {
  
  # Sanity checks
  if (!all(likelihood %in% c("gaussian", "bernoulli", "poisson")))
    stop("Liklihood not implemented: Use either gaussian, bernoulli or poisson")
  
  if (length(likelihood)==1) likelihood <- rep(likelihood, n_views) 
  if (!length(likelihood) == n_views) 
    stop("Likelihood needs to be a single string or matching the number of views!")
  
  # simulate facors
  Z <- matrix(rnorm(n_factors*n_samples, 0, 1), nrow=n_samples, ncol=n_factors)
  
  # set sparsity
  theta <- 0.5
  
  # set ARD prior, each factor being active in at least one view
  alpha <- vapply(seq_len(n_factors), function(fc) {
    active_vw <- sample(seq_len(n_views), 1)
    alpha_fc <- sample(c(1, 1000), n_views, replace = TRUE)
    if(all(alpha_fc==1000)) alpha_fc[active_vw] <- 1
    alpha_fc
  }, numeric(n_views))
  alpha <- matrix(alpha, nrow=n_factors, ncol=n_views, byrow=TRUE)
  
  # simulate weights 
  S <- lapply(seq_len(n_views), function(vw) matrix(rbinom(n_features*n_factors, 1, theta),
                                             nrow=n_features, ncol=n_factors))
  W <- lapply(seq_len(n_views), function(vw) vapply(seq_len(n_factors), function(fc) rnorm(n_features, 0, sqrt(1/alpha[fc,vw])), numeric(n_features)))
  
  # set noise level (for gaussian likelihood)
  tau <- 10
  
  # pre-compute linear term 
  mu <- lapply(seq_len(n_views), function(vw) Z %*% t(S[[vw]]*W[[vw]]))
  
  # simulate data according to the likelihood
  data <- lapply(seq_len(n_views), function(vw){
    lk <- likelihood[vw]
    if (lk == "gaussian"){
      dd <- t(mu[[vw]] + rnorm(length(mu[[vw]]),0,sqrt(1/tau)))
    }
    else if (lk == "poisson"){
      term <- log(1+exp(mu[[vw]]))
      dd <- t(apply(term, 2, function(tt) rpois(length(tt),tt)))
    }
    else if (lk == "bernoulli") {
      term <- 1/(1+exp(-mu[[vw]]))
      dd <- t(apply(term, 2, function(tt) rbinom(length(tt),1,tt)))
    }
    colnames(dd) <- paste0("sample_", seq_len(ncol(dd)))
    rownames(dd) <- paste0("feature", seq_len(nrow(dd)))
    dd
  })
  names(data) <- paste0("view_", seq_len(n_views))
  
  return(data)
}
