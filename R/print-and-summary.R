#' Print method for hmm objects
#' 
#' The \code{print} method for hmm objects displays a compact summary of the
#' fitted model. See the \strong{Details} section below for descriptions of the
#' different components of the printed output. For additional summary statistics
#' and diagnostics use the \code{\link[=summary.hmm]{summary}} method.
#' 
#' @export
#' @method print hmm 
#' @param digits Number of digits to use for formatting numbers.
#' @param ... Ignored.
#' @return Returns \code{x}, invisibly.  #' @details 
#' \subsection{Point estimates}{
#' Point estimates are medians computed from simulations.
#' For models fit using MCMC (\code{"sampling"}) the posterior
#' sample is used.  The point estimates reported are the same as the values
#' returned by \code{\link[=coef.stapreg]{coef}}.
#' }
#' \subsection{Uncertainty estimates (MAD_SD)}{
#' The standard deviations reported (labeled \code{MAD_SD} in the print output)
#' are computed from the same set of draws described above and are proportional
#' to the median absolute deviation (\code{\link[stats]{mad}}) from the median.
#' Compared to the raw posterior standard deviation, the MAD_SD will be
#' more robust for long-tailed distributions. 
#' }
#' 
#' @seealso \code{\link{summary.hmm}}, \code{\link{stapreg-methods}}
#' 
print.hmm <- function(x, digits = 2, ...) {
  cat("\n observations:", x$n)
  cat("\n groups:", x$J)

  cat("\n------\n")
  cat("Cluster Statistics")
  cat("\n------\n")

  pi_stats <- rbind(Median = apply(coda:::as.matrix.mcmc.list(x$pi),2,median),
                    MAD = apply(coda:::as.matrix.mcmc.list(x$pi),2,mad))

  .printfr(pi_stats,digits)

  mat <- cbind(min = min(sapply(x$num_clusters,min)),
            median = median(sapply(x$num_clusters,median)),
            max = max(sapply(x$num_clusters,max)))
  rownames(mat) <- "# of Clusters"

  .printfr(mat,digits)
  cat("\n------\n")

  w_stats <- rbind(Median = apply(coda:::as.matrix.mcmc.list(x$w),2,median),
                    MAD = apply(coda:::as.matrix.mcmc.list(x$w),2,mad))
  mu_stats <- rbind(Median = apply(coda:::as.matrix.mcmc.list(x$mu),2,median),
                    MAD = apply(coda:::as.matrix.mcmc.list(x$mu),2,mad))
  tau_stats <- rbind(Median = apply(coda:::as.matrix.mcmc.list(x$tau),2,median),
                     MAD = apply(coda:::as.matrix.mcmc.list(x$tau),2,mad))

  tau_switch <- all(dim(tau_stats) == dim(mu_stats))

  K <- max(x$if_df$Intensity_Function)
  L <- ncol(w_stats)/K

  for(k in 1:K){
      start <- (1 + L*(k-1))
      end <- (L + L*(k-1))
      if(pi_stats["Median",k]>.01){
          .printfr(w_stats[,start:end],digits)
          .printfr(mu_stats[,start:end],digits)
          if(tau_switch)
              .printfr(tau_stats[,start:end],digits)
      }
      cat("\n")
  }
  if(!tau_switch)
      .printfr(tau_stats,digits)

  cat("--- \n") 

  
  cat("* For help interpreting the printed output see ?print.hmm\n")
  
  invisible(x)
}


#' Summary method for hmm objects
#' 
#' Summaries of parameter estimates and MCMC convergence diagnostics 
#' (Monte Carlo error, effective sample size, Rhat).
#' 
#' @export
#' @method summary hmm 
#' 
#' 
#' @param ... Currently ignored.
#'   
#' @param probs For models fit using MCMC, 
#'   an optional numeric vector of probabilities passed to 
#'   \code{\link[stats]{quantile}}.
#' @param digits Number of digits to use for formatting numbers when printing. 
#'   When calling \code{summary}, the value of digits is stored as the 
#'   \code{"print.digits"} attribute of the returned object.
#'   
#' @return The \code{summary} method returns an object of class 
#'   \code{"summary.hmm"}.
#' 
#' 
#' @importFrom coda effectiveSize 
summary.hmm <- function(object,probs,digits = 1, ...) {


    summ <- function(x,y){ 
        out <- t(apply(coda:::as.matrix.mcmc.list(x[[y]]),2,function(a) c("Mean"=mean(a),"SD" = sd(a), quantile(a),"Geweke" = unname(coda::geweke.diag(a)$z) )))
        out <- cbind(out,"ESS"=effectiveSize(x[[y]]))
        if(length(x[[y]])>1)
            out <- cbind(out,"Rhat" = coda::gelman.diag(x[[y]],multivariate = FALSE,autoburnin=FALSE)$psrf[,1] )
        return(out)
        }
    out <- rbind(out,summ(object,"pi"),summ(object,"w"),summ(object,"mu"),summ(object,"tau"))


  structure(
    out,
    nobs = object$n,
    groups = object$J,
    posterior_sample_size = nrow(coda:::as.matrix.mcmc.list(object$alpha)),
    call = object$call,
    print.digits = digits,
    class = "summary.hmm"
  )

}


#' @rdname summary.hmm
#' @export
#' @method print summary.hmm
#'
#' @param x An object of class \code{"summary.hmm"}.
print.summary.hmm <- function(x, digits = max(1, attr(x, "print.digits")), 
                                  ...) {
    atts <- attributes(x)
    cat("\nModel Info:\n")
  cat("\n sample:      ", atts$posterior_sample_size, "(posterior sample size)")
  cat("\n observations:", atts$nobs)
  cat("\n groups:", atts$groups)

  cat("\n\nEstimates:\n")
    .printfr(x,digits)

    invisible(x)
}

#' converts summary output to latex
#' @export
#' @method to_latex hmm
#' @param object an hmm object
#' @param digits number of digits to round to
#' @param caption self-explanatory
#' 
to_latex <- function(object, digits = 1, caption = "")
	UseMethod("to_latex")

#' converts summary output to latex table
#' @export
#' 
to_latex.hmm <- function(object, digits = 1, caption ="", par = c("alpha")){


	s <- xtable::xtable(summary(object[par])$statistic,caption = caption, digits = digits)


    return(s)
}



# internal ----------------------------------------------------------------
.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}
