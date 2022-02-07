### Original genomescope model
nls_4peak <- function(x, y, k, estKmercov, estLength, max_iterations = 40){
    model4 = NULL

    cat("Original genomescope model 'nls_4peak'\n")

    try(model4 <- nls(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                          (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length +
                          (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length +
                          (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length),
                      start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = F)

    if(class(model4) == "try-error"){

        try(model4 <- nls(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                              (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length +
                              (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length +
                              (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length),
                          start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                          algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = F)
    }

    return(model4)
}

# two tissue model
nlsLM_2peak_unconditional_peaks <- function(x, y, kmerEst, lengthEst, hetEst = 0.6){
  nlsLM(y ~ ((het       * dnbinom(x, size = kmercov   / bias, mu = kmercov)) +
            ((1 - het)  * dnbinom(x, size = kmercov2  / bias, mu = kmercov2))) * length,
        start = list(kmercov = kmerEst, kmercov2 = (2 * kmerEst), bias = 0.5, length = lengthEst, het = hetEst),
        control = list(minFactor=1e-12, maxiter=40))
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # #  CONSIDERED MODELS # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # for now I will just try to create a model that will fit 1n and 2n peaks indipendently
# nls_2peak_kmer_explicit <- function(x, y, k, estKmercov, estLength, max_iterations = 40){
#     model2 = NULL
#
#   cat("trying nls_2peak_kmer_explicit standard algorithm\n")
#
#     try(model2 <- nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
#                           ((1-r)^k)        * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2)) * length,
#                       start = list(r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
#                       control = list(minFactor=1e-12, maxiter=max_iterations)), silent = F)
#
#     return(model2)
# }

# # for now I will just try to create a model that will fit 1n and 2n peaks indipendently
# nls_2peak_conditional <- function(x, y, k, estKmercov, estLength, rEst, biasEst, max_iterations = 40){
#     model2 = NULL
#
#     cat("trying nls_2peak_conditional standard algorithm\n")
#
#     try(model2 <- nlsLM(y ~ ((2*(1-(1-rEst)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
#                             ((1-rEst)^k)        * dnbinom(x, size = kmercov2  / bias, mu = kmercov2)) * estLength,
#                       start = list(kmercov=estKmercov, kmercov2 = (estKmercov * 2), bias = biasEst),
#                       control = list(minFactor=1e-12, maxiter=max_iterations)), silent = F)
#
#     return(model2)
# }
#
# # this model does not work well, for some reason it est negative fraction, I will try to remove the X chromosome term
# nls_2peak_two_tissue_constant_shift <- function(x, y, k, estKmercov, estLength, estR, biasEst, max_iterations = 40){
#     model2 = NULL
#
#     cat("trying nls_2peak_two_tissue_constant_shift standard algorithm\n")
#     # fixed parameters
#     r = estR
#     # length = estLength
#
#     try(model2 <- nlsLM(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov + cov_shift) +
#                             ((1-r)^k)        * dnbinom(x, size = (2 * kmercov) / bias, mu = (2 * kmercov) + cov_shift)) * length * fraction_diploid +
#                                                dnbinom(x, size = kmercov / bias, mu = kmercov + cov_shift) * length * (1 - fraction_diploid),
#                       start = list(kmercov=estKmercov, cov_shift=0, bias = biasEst, length = estLength, fraction_diploid=0.8),
#                       lower = c("kmercov" = 0, "cov_shift" = 0, "bias" = 0, "length" = 0, "fraction_diploid" = 0),
#                       upper = c("kmercov" = Inf, "cov_shift" = Inf, "bias" = Inf, "length" = Inf,  "fraction_diploid" = 1),
#                       control = list(minFactor=1e-12, maxiter=max_iterations)), silent = F)
#
#     return(model2)
# }
#
# # for now I will just try to create a model that will fit 1n and 2n peaks indipendently
# nls_4peak_two_tissue_constant_shift <- function(x, y, k, estKmercov, estLength, rEst, max_iterations = 40){
#     model4 = NULL
#
#     cat("trying nls_4peak_two_tissue_constant_shift standard algorithm\n")
#     # fixed paramters
#     r = rEst
#     # length = estLength
#     # fraction_diploid = fraction_diploid_est
#     # k as well of course
#
#     try(model4 <- nlsLM(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k)))  * dnbinom(x, size = (kmercov + cov_shift)   / bias, mu = kmercov + cov_shift)       +
#                              (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = (kmercov + cov_shift)*2 / bias, mu = (kmercov * 2) + cov_shift) +
#                              (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = (kmercov + cov_shift)*3 / bias, mu = (kmercov * 3) + cov_shift) +
#                              (d*(1-r)^(2*k))                                                              * dnbinom(x, size = (kmercov + cov_shift)*4 / bias, mu = (kmercov * 4) + cov_shift)) * length * fraction_diploid +
#                             ((1 - d)                                                                      * dnbinom(x, size = (kmercov + cov_shift)   / bias, mu = kmercov + cov_shift)       +
#                              d                                                                            * dnbinom(x, size = (kmercov + cov_shift)*2 / bias, mu = (kmercov * 2) + cov_shift)) * length * (1 - fraction_diploid),
#                         start = list(d=0, kmercov=estKmercov, bias = 0.5, cov_shift = 0, fraction_diploid = 0.7, length = estLength),
#                         control = list(minFactor=1e-12, maxiter=max_iterations)), silent = F)
#
#     return(model4)
# }
