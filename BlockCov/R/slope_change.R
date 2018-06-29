#' This function fits to a numerical vector sorted in the non decreasing order two simple linear regressions and returns the index corresponding to the estimated change between the two regression models.
#'
#' @param  Y numerical vector sorted in the non decreasing order.
#' @return K the index corresponding to the estimated change between the two linear regression models.
#' @importFrom Matrix Matrix
#' @examples
#' n <- 30
#' q <- 100
#' Sigma <- Simu_Sigma(q = q, diag = FALSE, equal = TRUE)
#' Matrix::image(Sigma)
#' E <- matrix(rnorm(n * q), ncol = q) %*% chol(as.matrix(Sigma))
#' corE <- cor(as.matrix(E))
#' vec_up_emp <- corE[upper.tri(corE)]
#' Pti_sig <- matrix(0, ncol = (q - 1), nrow = (q - 1))
#' Pti_sig[upper.tri(Pti_sig, diag=TRUE)] <- vec_up_emp
#' Pti_sig[lower.tri(Pti_sig)] <- t(as.matrix(Pti_sig))[lower.tri(t(as.matrix(Pti_sig)))]
#' res_svd <- svd(Pti_sig)
#' vp      <- res_svd$d
#' slope_change(vp)
#' @export
slope_change <- function(Y){
  x=1:length(Y);
  l=length(Y);
  kvect=2:(l-1);

  achapL=rep(0,length(kvect));
  bchapL=rep(0,length(kvect));
  achapR=rep(0,length(kvect));
  bchapR=rep(0,length(kvect));
  S=rep(0,length(kvect));

  for (j in 1:length(kvect)){
    k=kvect[j];
    YL=Y[1:k];
    xL=1:k;

    achapL[k]=(mean(YL*xL)-mean(YL)*mean(xL))/(mean(xL^2)-(mean(xL))^2);
    bchapL[k]=mean(YL)-achapL[k]*mean(xL);

    YR=Y[k:length(Y)];
    xR=k:length(Y);

    achapR[k]=(mean(YR*xR)-mean(YR)*mean(xR))/(mean(xR^2)-(mean(xR))^2);
    bchapR[k]=mean(YR)-achapR[k]*mean(xR);

    S[j]=sum((YL-achapL[k]*xL-bchapL[k])^2)+sum((YR-achapR[k]*xR-bchapR[k])^2);
  }
  K=kvect[which.min(S)];
  return(K);
}
