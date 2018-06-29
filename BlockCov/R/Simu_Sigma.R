#' This function generates a block structured symmetric positive definite matrix to test the BlockCov methodology.
#'
#' @param  q integer corresponding to the size of the covariance matrix.
#' @param  diag logical, whether or not the covariance matrix is block-diagonal.
#' @param  equal logical, whether or not the values in the blocks are equal.
#' @return Sigma a correlation matrix to test the BlockCov methodology.
#' @importFrom Matrix Matrix
#' @examples
#'Sigma <- Simu_Sigma(q = 100, diag = FALSE, equal = TRUE)
#'Matrix::image(Sigma)
#' @export
Simu_Sigma <- function(q, diag = TRUE, equal = TRUE){

  list_a=c(floor(0.1*q),floor(0.2*q),floor(0.3*q),floor(0.2*q),floor(0.2*q))
  list_rho=c(0.7,0.75,0.65,0.8,0.7)

  nb_bloc=length(list_a)
  position=cumsum(list_a)
  Z=Matrix(0,nrow=q,ncol=nb_bloc)
  if (equal){
    Z[1:position[1],1]=rep(sqrt(list_rho[1]),list_a[1])
    for (i in 2:nb_bloc){
      Z[(position[(i-1)]+1):position[i], i] = rep(sqrt(list_rho[i]), list_a[i])
    }
  }else{
    Z[1:position[1],1]=runif(list_a[1],sqrt(0.6),sqrt(0.8))
    for (i in 2:nb_bloc){
      if(i == 3) {
        Z[(position[(i-1)]+1):position[i], i] = runif(list_a[i],sqrt(0.3),sqrt(0.4))
      } else{
        Z[(position[(i-1)]+1):position[i], i] = runif(list_a[i],sqrt(0.6),sqrt(0.8))
      }  }
  }

  if (!diag) Z[floor(0.35*q):floor(0.45*q),4]=-0.5
  Sigma=Z%*%t(Z)
  diag(Sigma)=1
  return(Sigma)
}

