
# Print the results.
#' @export
print.trm <- function(x, ...){
  cat("Triadic Relations Model, analysized by triadR (Wong & Kenny, 2018)\n")
  cat("------------------------------------------------------------------------\n")
  if (x$is.univar == TRUE){
    var.results <- x$varComp[1:7,-5]
    colnames(var.results)[2] <- "proportion"
    cov.results <- x$varComp[8:23,-5]
    colnames(cov.results)[2] <- "correlation"
    cat("Number of group    :  ",nrow(x$varComp.raw),"\n")
    cat("Average group size :  ", mean(x$varComp.raw$n),"\n")
    if (nrow(x$varComp.raw)<2){
      x$varComp <- x$varComp[,c(1,2)]
      cat("\nNote:\tSignificance tests are not available for single group analysis.\n")
    }
    if (any(abs(cov.results[,2]) > 1)){
      cat("\nNote:\tIt is possible for correlations computed using method-of-moments\n")
      cat("\tto produce values larger than one in absolute value. \n")
      cat("\tIf such values occur, they would typically be reported as either\n")
      cat("\t1.00 or -1.00.\n")
    }
    cat("------------------------------------------------------------------------\n")
    cat("Univariate analysis for  : ", x$var.1, "\n\n")
    cat("Variance Component:\n")
    print(round(var.results, digits = 3))
    cat("\nCovariance Component:\n")
    cat("\nIndividual\n")
    print(round(cov.results[1:3,], digits = 3))
    cat("\nDyadic\n")
    print(round(cov.results[4:12,], digits = 3))
    cat("\nTriadic\n")
    print(round(cov.results[13:16,], digits = 3))
  }else{
    cat("Number of group    :  ",nrow(x$varComp.raw$var.1),"\n")
    cat("Average group size :  ", mean(x$varComp.raw$var.1$n),"\n")
    if (nrow(x$varComp.raw$var.1)<2){
      x$varComp$var.1 <- x$varComp$var.1[,c(1,2)]
      x$varComp$var.2 <- x$varComp$var.2[,c(1,2)]
      x$varComp$bivar <- x$varComp$bivar[,c(1,2)]
      cat("\nNote:\tSignificance tests are not available for single group analysis.\n")
    }
    cat("\nNote:\tIt is possible for correlations computed using method-of-moments\n")
    cat("\tto produce values larger than one in absolute value. \n")
    cat("\tIf such values occur, they would typically be reported as either\n")
    cat("\t1.00 or -1.00.\n")
    cat("------------------------------------------------------------------------\n")
    cat("Univariate analysis for  : ", x$var.1, "\n\n")
    cat("Variance Component:\n")
    print(round(x$varComp$var.1[1:7,], digits = 3))
    cat("\nCovariance Component:\n")
    cat("\nIndividual\n")
    print(round(x$varComp$var.1[8:10,], digits = 3))
    cat("\nDyadic\n")
    print(round(x$varComp$var.1[11:19,], digits = 3))
    cat("\nTriadic\n")
    print(round(x$varComp$var.1[20:23,], digits = 3))
    cat("------------------------------------------------------------------------\n")
    cat("Univariate analysis for  : ", x$var.2, "\n\n")
    cat("Variance Component:\n")
    print(round(x$varComp$var.2[1:7,], digits = 3))
    cat("\nCovariance Component:\n")
    cat("\nIndividual\n")
    print(round(x$varComp$var.2[8:10,], digits = 3))
    cat("\nDyadic\n")
    print(round(x$varComp$var.2[11:19,], digits = 3))
    cat("\nTriadic\n")
    print(round(x$varComp$var.2[20:23,], digits = 3))
    cat("------------------------------------------------------------------------\n")
    cat("Bivariate analysis for   : ", x$var.1, " and ", x$var.2, "\n\n")
    cat("Covariance Component:\n")
    cat("\nIndividual\n")
    print(round(x$varComp$bivar[1:9,], digits = 3))
    cat("\nDyadic\n")
    print(round(x$varComp$bivar[10:27,], digits = 3))
    cat("\nTriadic\n")
    print(round(x$varComp$bivar[28:33,], digits = 3))
  }

}
