### Original genomescope model
nls_4peak <- function(x, y, k, estKmercov, estLength, max_iterations){
    model4 = NULL

    cat("Fitting patemeters: \n")
    cat(paste("Fitting coverage range: ", min(x), max(x), '\n'))
    cat(paste("Fitting frequency range: ", min(y), max(y), '\n'))
    cat(paste0("\tk:\t", k, "\n"))
    cat(paste0("\testKmercov:\t", estKmercov, "\n"))
    cat(paste0("\testLength:\t", estLength, "\n"))
    cat(paste0("\tmax_iterations:\t", max_iterations, "\n"))

    try(model4 <- nls(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                          (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length +
                          (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length +
                          (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length),
                      start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    if(class(model4) == "try-error"){

        try(model4 <- nls(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                              (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length +
                              (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length +
                              (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length),
                          start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                          algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
    }

    return(model4)
}



# for now I will just try to create a model that will fit 1n and 2n peaks indipendently

nls_2peak_kmer_explicit <- function(x, y, k, estKmercov, estLength, max_iterations){
    model2 = NULL

  cat("trying nls_2peak_kmer_explicit standard algorithm\n")

    try(model2 <- nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                        ((1-r)^k)        * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2)) * length,
                      start = list(r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}

# for now I will just try to create a model that will fit 1n and 2n peaks indipendently

nls_2peak_conditional <- function(x, y, k, estKmercov, estLength, rEst, biasEst, max_iterations){
    model2 = NULL

    cat("trying nls_2peak_conditional standard algorithm\n")

    try(model2 <- nls(y ~ ((2*(1-(1-rEst)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                          ((1-rEst)^k)        * dnbinom(x, size = kmercov2  / bias, mu = kmercov2)) * estLength,
                      start = list(kmercov=estKmercov, kmercov2 = (estKmercov * 2), bias = biasEst),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}


nls_2peak_X0_conditional <- function(x, y, k, estKmercov, estLength, rEst, biasEst, max_iterations){
    model2 = NULL

    cat("trying nls_2peak_conditional standard algorithm\n")

    try(model2 <- nls(y ~ ((2*(1-(1-rEst)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                          ((1-rEst)^k)        * dnbinom(x, size = kmercov2  / bias, mu = kmercov2)) * estLength,
                      start = list(kmercov=estKmercov, kmercov2 = (estKmercov * 2), bias = biasEst),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}

nls_2peak_male_X0_conditional <- function(x, y, k, estKmercov, estLength, estR, biasEst, max_iterations){
    model2 = NULL

    cat("trying nls_2peak_male_X0_conditional standard algorithm\n")
    # fixed parameters
    r = estR
    length = estLength

    try(model2 <- nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                          ((1-r)^k)        * dnbinom(x, size = kmercov2 / bias, mu = kmercov2)) * length * fraction_diploid +
                          dnbinom(x, size = kmercov / bias, mu = kmercov) * length * (1 - fraction_diploid),
                      start = list(kmercov=estKmercov, kmercov2=(estKmercov * 2), bias = biasEst, fraction_diploid=0.8),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}


nls_2peak_two_tissue_constant_shift <- function(x, y, k, estKmercov, estLength, estR, biasEst, max_iterations){
    model2 = NULL

    cat("trying nls_2peak_two_tissue_constant_shift standard algorithm\n")
    # fixed parameters
    r = estR
    length = estLength

    try(model2 <- nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov + cov_shift) +
                          ((1-r)^k)        * dnbinom(x, size = (2 * kmercov) / bias, mu = (2 * kmercov) + cov_shift)) * length * fraction_diploid +
                          dnbinom(x, size = kmercov / bias, mu = kmercov + cov_shift) * length * (1 - fraction_diploid),
                      start = list(kmercov=estKmercov, cov_shift=0, bias = biasEst, fraction_diploid=0.8),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}

# PGE model
nls_3peak_two_tissue_PGE <- function(x, y, k, estKmercov, estLength, estR, biasEst, estFrac, max_iterations){
    model2 = NULL

    cat("trying nls_2peak_two_tissue_constant_shift standard algorithm\n")
    # fixed parameters
    r = estR                 # autosomal heterozygosity
    # length = estLength       # haploid genome length
    # fraction_diploid = estFrac # estimated fraction of Autosomes

    try(model2 <- nls(y ~ ((1-(1-r)^k) * dnbinom(x, size = kmercov_p / bias, mu = kmercov_p) +                                                             # paternal heterozygous autosomes
                           (1-(1-r)^k) * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m) +                                                             # maternal heterozygous autosomes
                           ((1-r)^k)   * dnbinom(x, size = (kmercov_p + kmercov_m) / bias, mu = kmercov_p + kmercov_m)) * length * fraction_diploid +      # homozygous autosomes
                           1           * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m)                            * length * (1 - fraction_diploid), # maternally inherited X autosomes
                      start = list(kmercov_p=estKmercov, kmercov_m=estKmercov, bias=biasEst, length=estLength, fraction_diploid = estFrac),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = F)

    return(model2)
}

# for now I will just try to create a model that will fit 1n and 2n peaks indipendently
# I changed just tiny buits

nls_2peak_general <- function(x, y, estKmercov, max_iterations){

    model2 <- NULL

  cat("trying nls_2peak_general standard algorithm\n")
    try(model2 <- nls(y ~ N * (dnbinom(x, size = kmercov_1n / bias, mu = kmercov_1n) +
                               dnbinom(x, size = kmercov_2n / bias, mu = kmercov_2n)),
                    start = list(kmercov_1n=estKmercov, kmercov_2n=(2 * estKmercov), bias = 0.5, N = sum(y)),
                    control = list(minFactor=1e-12, maxiter=max_iterations)))

    if(is.null(model2)){
        cat("retrying nls_2peak_general with port algorithm\n")
        try(model2 <- nls(y ~ N * (dnbinom(x, size = kmercov_1n / bias1, mu = kmercov_1n) +
                                         dnbinom(x, size = kmercov_2n / bias2, mu = kmercov_2n)),
                                start = list(kmercov_1n=estKmercov, kmercov_2n=(2 * estKmercov), bias1 = 1, bias2 = 2, N = sum(y)),
                          algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)))
    }

    return(model2)
}


predict_maternal_3peak_PGE <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  ((1-(1-r)^k)     * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m) +
  ((1-r)^k) * 1/2  * dnbinom(x, size = (kmercov_p + kmercov_m) / bias, mu = kmercov_p + kmercov_m)) * length * fraction_diploid +
  1                * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m)                            * length * (1 - fraction_diploid)
}

predict_maternal_3peak_PGE_het <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  ((1-(1-r)^k)     * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m)) * length * fraction_diploid
}

predict_maternal_3peak_PGE_autosomes <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  ((1-(1-r)^k)     * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m) +
  ((1-r)^k) * 1/2  * dnbinom(x, size = (kmercov_p + kmercov_m) / bias, mu = kmercov_p + kmercov_m)) * length * fraction_diploid
}

predict_maternal_3peak_PGE_X <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
   dnbinom(x, size = kmercov_m / bias, mu = kmercov_m) * length * (1 - fraction_diploid)
}

predict_paternal_3peak_PGE <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  ((1-(1-r)^k)     * dnbinom(x, size = kmercov_p / bias, mu = kmercov_p) +
  ((1-r)^k) * 1/2 * dnbinom(x, size = (kmercov_p + kmercov_m) / bias, mu = kmercov_p + kmercov_m)) * length * fraction_diploid
}

predict_paternal_3peak_PGE_het <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  ((1-(1-r)^k)     * dnbinom(x, size = kmercov_p / bias, mu = kmercov_p)) * length * fraction_diploid
}

predict_shared_3peak_PGE <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  (((1-r)^k) * dnbinom(x, size = (kmercov_p + kmercov_m) / bias, mu = kmercov_p + kmercov_m)) * length * fraction_diploid
}

plot_nls_3peak_two_tissue_PGE_model(x, y, est_model, estR){

  fitted_parameters <- coef(est_model)
  cov_m <- fitted_parameters['kmercov_m']
  cov_p <- fitted_parameters['kmercov_p']
  cov <- cov_m + cov_p

  matermal_prediction <- predict_maternal_3peak_PGE(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  matermal_prediction_het <- predict_maternal_3peak_PGE_het(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  matermal_prediction_A <- predict_maternal_3peak_PGE_autosomes(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  matermal_prediction_X <- predict_maternal_3peak_PGE_X(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  paternal_prediction <- predict_paternal_3peak_PGE(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  paternal_prediction_het <- predict_paternal_3peak_PGE_het(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  shared_prediction <- predict_shared_3peak_PGE(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])

  barplot <- barplot(y ~ x, col = 'deepskyblue', border = F, xlab = 'Coverage', ylab = 'Frequency', ylim = c(0, max(y) * 1.1))

  lines(predict(est_model, response = T) ~ barplot, lwd = 3)

  # lines(matermal_prediction ~ barplot, lwd = 3, col = 'darkorange2')
  lines(matermal_prediction_het ~ barplot, lwd = 3, col = 'firebrick3')
  lines(paternal_prediction_het ~ barplot, lwd = 3, col = 'darkblue')
  lines(shared_prediction ~ barplot, lwd = 3, col = 'darkorchid3')
  lines(matermal_prediction_X ~ barplot, lwd = 3, col = "gold")


  # lines(c(cov_m, cov_m), c(-1e10, 1e10), lwd = 2, lty = 2, col = 'firebrick3')
  # lines(c(cov_p, cov_p), c(-1e10, 1e10), lwd = 2, lty = 2, col = 'darkblue')
  # lines(c(cov, cov), c(-1e10, 1e10), lwd = 2, lty = 2, col = 'darkorchid3')

  legend('topright',
         c('kmer histogram','full model', 'shared A', 'maternal X chromosomes', 'maternal heterozygous A', 'paternal heterozygous A'),
         col = c('deepskyblue','black', 'darkorchid3', 'gold', 'firebrick3', 'darkblue'),
         lty = c(NA, 1, 1, 1, 1, 1), lwd = 3,
         pch = c(15, NA, NA, NA, NA, NA), bty = 'n')

  title(paste0('Explicitly fitted PGE model to A. fusca'))

}
