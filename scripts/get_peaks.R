get_peaks <- function(ks){
  second_deriv <- diff(sign(diff(ks$y)))

  peak_covs <- ks$x[which(second_deriv == -2) + 1]
  peak_heights <- ks$y[which(second_deriv == -2) + 1]
  peaks <- data.frame(cov = peak_covs, height = peak_heights)
	peaks[order(peaks$height), ]
}
