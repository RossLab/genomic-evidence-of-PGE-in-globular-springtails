require('minpack.lm')

###############
## FUNCTIONS ##
###############

### Original genomescope model
nls_4peak <- function(x, y, k, estKmercov, estLength, max_iterations = 40){
    model4 = NULL

    cat("Original genomescope model 'nls_4peak'\n")

    try(model4 <- nlsLM(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
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

predict_1n_peak_nlsLM_2peak_unconditional <- function(model){
  model_env <- model$m$getEnv()
  model_env$het * dnbinom(model_env$x, size = model_env$kmercov   / model_env$bias, mu = model_env$kmercov) * model_env$length
}

predict_2n_peak_nlsLM_2peak_unconditional <- function(model){
  model_env <- model$m$getEnv()
  (1 - model_env$het)  * dnbinom(model_env$x, size = model_env$kmercov2  / model_env$bias, mu = model_env$kmercov2) * model_env$length
}

coverage_barplot <- function(bar_heights, bar_positions, font_size = 1, width = 0.5){

  plot(bar_heights, type="n", xlab="Coverage", ylab="Frequency",
       ylim=c(0, max(bar_heights)), xlim=range(bar_positions),
       cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  for ( i in 1:length(bar_heights)){
    rect(bar_positions[i] - width, 0, bar_positions[i] + width, bar_heights[i], col = 'deepskyblue', border = F)
  }
}

plot_PGE_model <- function(model){
  x <- model$m$getEnv()$x
  y <- model$m$getEnv()$y
  cov_1n <- coef(model)['kmercov']
  cov_2n <- coef(model)['kmercov2']
  monoploid_col <- 'darkorchid4'
  diploid_col <- 'chocolate'

  coverage_barplot(y, x)

  disomic_prediction <- predict_2n_peak_nlsLM_2peak_unconditional(model)
  monosomic_prediction <- predict_1n_peak_nlsLM_2peak_unconditional(model)

  lines(predict(model, response = T) ~ x, lwd = 3)
  lines(monosomic_prediction ~ x, lwd = 3, col = monoploid_col)
  lines(disomic_prediction ~ x, lwd = 3, col = diploid_col)

  lines(c(cov_1n, cov_1n) - 0.5, c(0, max(y)), lty = 2)
  lines(c(cov_2n, cov_2n) - 0.5, c(0, max(y)), lty = 2)

  legend('topright',
         c('kmer histogram','full model', 'monoploid', 'diploid'),
         col = c('deepskyblue','black', monoploid_col, diploid_col),
         lty = c(NA, 1, 1, 1), lwd = 3,
         pch = c(15, NA, NA, NA), bty = 'n')
}
