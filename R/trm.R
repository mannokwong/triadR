#' triadR: Analyzing Generalized Round-Robin Data Using R
#'
#' The function trm() performs Tradic Relations Model (TRM) analyses for single or multiple generalized round-robin groups.
#'
#' @param formula formula
#' @param data data.frame in long format
#' @examples
#' # Load the package
#' library(triadR)
#'
#' # Univariate Analysis
#' data(PM1997)
#' trm.1 <- trm(var ~ jid*aid*pid, PM1997)
#' @export

trm <- function(formula, data) {
  # Read the formula
  sf  <- trmGlobal$suffixes
  FC <- trm.check(formula)

  # Check the dataset;
  # Create generalized RR matrices
  GRR <- trm.create.grr(FC, data)

  if(FC[[7]]==TRUE){ ### Univariate Analysis
    var1.grr.means  <- trm.grr.means(GRR[[1]])
    var1.effects  <- trm.effects(var1.grr.means, FC[[5]], FC[[4]], FC[[1]], FC[[2]], FC[[3]])

    rm(GRR)
    var1.cor  <- trm.univariate.cor(var1.grr.means, FC[[5]])
    var1.pool <- trm.pool(var1.cor, 1)

    names(var1.effects) <- c(paste(FC[[5]],sf[1], sep = ""),
                             paste(FC[[5]],sf[2], sep = ""),
                             paste(FC[[5]],sf[3], sep = ""),
                             paste(FC[[5]],sf[4], sep = ""),
                             paste(FC[[5]],sf[5], sep = ""),
                             paste(FC[[5]],sf[6], sep = ""),
                             paste(FC[[5]],sf[10], sep = ""))

    res  <- list(is.univar = FC[[7]], var.1 = FC[[5]], effects = var1.effects, varComp.raw = var1.cor, varComp = var1.pool[[1]])

  }else{ ### Bivariate Analysis
    var1.grr.means  <- trm.grr.means(GRR[[1]])
    var2.grr.means  <- trm.grr.means(GRR[[2]])
    var1.effects  <- trm.effects(var1.grr.means, FC[[5]], FC[[4]], FC[[1]], FC[[2]], FC[[3]])
    var2.effects  <- trm.effects(var2.grr.means, FC[[6]], FC[[4]], FC[[1]], FC[[2]], FC[[3]])
    rm(GRR)
    var1.cor  <- trm.univariate.cor(var1.grr.means, FC[[5]])
    var2.cor  <- trm.univariate.cor(var2.grr.means, FC[[6]])
    bivar.cor <- trm.bivariate.cor(var1.grr.means,var2.grr.means, FC[[5]], FC[[6]])

    var1.pool <- trm.pool(var1.cor, 1)
    var2.pool <- trm.pool(var2.cor, 1)
    bivar.pool <- trm.pool(bivar.cor, 2, var1.pool[[2]], var2.pool[[2]])

    names(var1.effects) <- c(paste(FC[[5]],sf[1], sep = ""),
                             paste(FC[[5]],sf[2], sep = ""),
                             paste(FC[[5]],sf[3], sep = ""),
                             paste(FC[[5]],sf[4], sep = ""),
                             paste(FC[[5]],sf[5], sep = ""),
                             paste(FC[[5]],sf[6], sep = ""),
                             paste(FC[[5]],sf[10], sep = ""))
    names(var2.effects) <- c(paste(FC[[6]],sf[1], sep = ""),
                             paste(FC[[6]],sf[2], sep = ""),
                             paste(FC[[6]],sf[3], sep = ""),
                             paste(FC[[6]],sf[4], sep = ""),
                             paste(FC[[6]],sf[5], sep = ""),
                             paste(FC[[6]],sf[6], sep = ""),
                             paste(FC[[6]],sf[10], sep = ""))

    res  <- list(is.univar = FC[[7]], var.1 = FC[[5]], var.2 = FC[[6]], effects = append(var1.effects, var2.effects), varComp.raw = list(var.1 = var1.cor, var.2 = var2.cor, bivar = bivar.cor), varComp = list(var.1 = var1.pool[[1]], var.2 = var2.pool[[1]], bivar = bivar.pool))

  }
  class(res) <- "trm"

  return(res)
}

trm.pool <- function(var.cor, idx, var1.estimates = NULL, var2.estimates = NULL){
  if (idx == 1){
    vars <- lapply(var.cor[,4:10],
                   function(x, n) {
                     temp <- data.frame(x,n)
                     temp <- temp[!is.na(temp[,1]),]
                     ngroups <- length(x[!is.na(x)])
                     #Assuming weights are n - 1
                     temp$weight <- temp$n-1
                     sum_weight <- sum(temp$weight)
                     #Relative weights
                     temp$rel_weight <- temp$weight/sum_weight

                     w_mean <- weighted.mean(temp$x,temp$rel_weight)
                     wvar   <- ngroups*sum(temp$rel_weight*(temp$x-w_mean)^2)/(ngroups-1)
                     ttest  <- w_mean*ngroups/sqrt(wvar)
                     df     <- ngroups - 1
                     # pvalue for a variance, one-tailed
                     pval   <- 1-pt(ttest, ngroups-1, 0, lower.tail = TRUE, log.p = FALSE)
                     res    <- c(w_mean, sqrt(wvar), ttest, ngroups, df, pval)
                     return(res)
                    }, var.cor[,2])


    covs <- lapply(var.cor[,11:26],
                   function(x, n) {
                     temp <- data.frame(x,n)
                     temp <- temp[!is.na(temp[,1]),]
                     ngroups <- length(x[!is.na(x)])
                     #Assuming weights are n - 1
                     temp$weight <- temp$n-1
                     sum_weight <- sum(temp$weight)
                     #Relative weights
                     temp$rel_weight <- temp$weight/sum_weight

                     w_mean <- weighted.mean(temp$x,temp$rel_weight)
                     wvar   <- ngroups*sum(temp$rel_weight*(temp$x-w_mean)^2)/(ngroups-1)
                     ttest  <- w_mean*ngroups/sqrt(wvar)
                     df     <- ngroups -1
                     # pvalue for a covariance, two-tailed
                     pval   <- 2*(1-pt(ttest, ngroups-1, 0, lower.tail = TRUE, log.p = FALSE))
                     res    <- c(w_mean, sqrt(wvar), ttest, ngroups, df, pval)
                     return(res)
                    }, var.cor[,2])
    var.estimates <- sapply(vars, function(x){return(x[1])})
    var.results   <- data.frame(matrix(unlist(vars), nrow=7, byrow=TRUE))
    var.results[,7] <- var.results[,1]/sum(var.estimates, na.rm = TRUE)

    cov.results   <- data.frame(matrix(unlist(covs), nrow=16, byrow=TRUE))
    cov.results[1,7] <- cov.results[1,1]/sqrt(var.estimates[1]*var.estimates[2])
    cov.results[2,7] <- cov.results[2,1]/sqrt(var.estimates[1]*var.estimates[3])
    cov.results[3,7] <- cov.results[3,1]/sqrt(var.estimates[2]*var.estimates[3])

    cov.results[4,7] <- cov.results[4,1]/sqrt(var.estimates[4]*var.estimates[4])
    cov.results[5,7] <- cov.results[5,1]/sqrt(var.estimates[5]*var.estimates[5])
    cov.results[6,7] <- cov.results[6,1]/sqrt(var.estimates[6]*var.estimates[6])

    cov.results[7,7] <- cov.results[7,1]/sqrt(var.estimates[4]*var.estimates[5])
    cov.results[8,7] <- cov.results[8,1]/sqrt(var.estimates[4]*var.estimates[6])
    cov.results[9,7] <- cov.results[9,1]/sqrt(var.estimates[5]*var.estimates[6])

    cov.results[10,7] <- cov.results[10,1]/sqrt(var.estimates[4]*var.estimates[5])
    cov.results[11,7] <- cov.results[11,1]/sqrt(var.estimates[4]*var.estimates[6])
    cov.results[12,7] <- cov.results[12,1]/sqrt(var.estimates[5]*var.estimates[6])

    cov.results[13,7] <- cov.results[13,1]/sqrt(var.estimates[7]^2)
    cov.results[14,7] <- cov.results[14,1]/sqrt(var.estimates[7]^2)
    cov.results[15,7] <- cov.results[15,1]/sqrt(var.estimates[7]^2)
    cov.results[16,7] <- cov.results[16,1]/sqrt(var.estimates[7]^2)

    res  <- rbind(var.results[,c(1,7,2,3,4,5,6)],cov.results[,c(1,7,2,3,4,5,6)])
    rownames(res) <- c(names(vars), names(covs))
    colnames(res) <- c("estimate", "standardized","se", "t.value", "n", "df", "p.value")
    return(list(res, var.estimates))

  }else if (idx == 2){
    covs <- lapply(var.cor[,3:35],
                   function(x, n) {
                     temp <- data.frame(x,n)
                     temp <- temp[!is.na(temp[,1]),]
                     ngroups <- length(x[!is.na(x)])
                     #Assuming weights are n - 1
                     temp$weight <- temp$n-1
                     sum_weight <- sum(temp$weight)
                     #Relative weights
                     temp$rel_weight <- temp$weight/sum_weight

                     w_mean <- weighted.mean(temp$x,temp$rel_weight)
                     wvar   <- ngroups*sum(temp$rel_weight*(temp$x-w_mean)^2)/(ngroups-1)
                     ttest  <- w_mean*ngroups/sqrt(wvar)
                     df     <- ngroups -1
                     # pvalue for a covariance, two-tailed
                     pval   <- 2*(1-pt(ttest, ngroups-1, 0, lower.tail = TRUE, log.p = FALSE))
                     res    <- c(w_mean, sqrt(wvar), ttest, ngroups, df, pval)
                     return(res)
                   }, var.cor[,2])
    cov.results   <- data.frame(matrix(unlist(covs), nrow=33, byrow=TRUE))

    cov.results[1,7] <- cov.results[1,1]/sqrt(var1.estimates[1]*var2.estimates[1])
    cov.results[2,7] <- cov.results[2,1]/sqrt(var1.estimates[2]*var2.estimates[2])
    cov.results[3,7] <- cov.results[3,1]/sqrt(var1.estimates[3]*var2.estimates[3])
    cov.results[4,7] <- cov.results[4,1]/sqrt(var1.estimates[1]*var2.estimates[2])
    cov.results[5,7] <- cov.results[5,1]/sqrt(var1.estimates[1]*var2.estimates[3])
    cov.results[6,7] <- cov.results[6,1]/sqrt(var1.estimates[2]*var2.estimates[1])
    cov.results[7,7] <- cov.results[7,1]/sqrt(var1.estimates[2]*var2.estimates[3])
    cov.results[8,7] <- cov.results[8,1]/sqrt(var1.estimates[3]*var2.estimates[1])
    cov.results[9,7] <- cov.results[9,1]/sqrt(var1.estimates[3]*var2.estimates[2])

    cov.results[10,7] <- cov.results[10,1]/sqrt(var1.estimates[4]*var2.estimates[4])
    cov.results[11,7] <- cov.results[11,1]/sqrt(var1.estimates[5]*var2.estimates[5])
    cov.results[12,7] <- cov.results[12,1]/sqrt(var1.estimates[6]*var2.estimates[6])
    cov.results[13,7] <- cov.results[13,1]/sqrt(var1.estimates[4]*var2.estimates[4])
    cov.results[14,7] <- cov.results[14,1]/sqrt(var1.estimates[5]*var2.estimates[5])
    cov.results[15,7] <- cov.results[15,1]/sqrt(var1.estimates[6]*var2.estimates[6])
    cov.results[16,7] <- cov.results[16,1]/sqrt(var1.estimates[4]*var2.estimates[5])
    cov.results[17,7] <- cov.results[17,1]/sqrt(var1.estimates[4]*var2.estimates[6])
    cov.results[18,7] <- cov.results[18,1]/sqrt(var1.estimates[5]*var2.estimates[7])
    cov.results[19,7] <- cov.results[19,1]/sqrt(var1.estimates[5]*var2.estimates[6])
    cov.results[20,7] <- cov.results[20,1]/sqrt(var1.estimates[6]*var2.estimates[4])
    cov.results[21,7] <- cov.results[21,1]/sqrt(var1.estimates[6]*var2.estimates[5])
    cov.results[22,7] <- cov.results[22,1]/sqrt(var1.estimates[4]*var2.estimates[5])
    cov.results[23,7] <- cov.results[23,1]/sqrt(var1.estimates[4]*var2.estimates[6])
    cov.results[24,7] <- cov.results[24,1]/sqrt(var1.estimates[5]*var2.estimates[4])
    cov.results[25,7] <- cov.results[25,1]/sqrt(var1.estimates[5]*var2.estimates[6])
    cov.results[26,7] <- cov.results[26,1]/sqrt(var1.estimates[6]*var2.estimates[4])
    cov.results[27,7] <- cov.results[27,1]/sqrt(var1.estimates[6]*var2.estimates[5])

    cov.results[28,7] <- cov.results[28,1]/sqrt(var1.estimates[7]*var2.estimates[7])
    cov.results[29,7] <- cov.results[29,1]/sqrt(var1.estimates[7]*var2.estimates[7])
    cov.results[30,7] <- cov.results[30,1]/sqrt(var1.estimates[7]*var2.estimates[7])
    cov.results[31,7] <- cov.results[31,1]/sqrt(var1.estimates[7]*var2.estimates[7])
    cov.results[32,7] <- cov.results[32,1]/sqrt(var1.estimates[7]*var2.estimates[7])
    cov.results[33,7] <- cov.results[33,1]/sqrt(var1.estimates[7]*var2.estimates[7])

    res  <- cov.results[,c(1,7,2,3,4,5,6)]
    rownames(res) <- names(covs)
    colnames(res) <- c("estimate", "standardized","se", "t.value", "n", "df", "p.value")
    return(res)

    ####
  }
}

trm.bivariate.cor <- function(var1.grr.means, var2.grr.means, var.1, var.2){

  cor.list  <- list()
  ### Calculate correlations for each group
  for (i in names(var1.grr.means)){
    temp.1  <- var1.grr.means[[i]]
    temp.2  <- var2.grr.means[[i]]

    cp.l_l  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[3]]-temp.1[[2]])*(temp.2[[3]]-temp.2[[2]]), na.rm = TRUE))
    cp.r_r  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[4]]-temp.1[[2]])*(temp.2[[4]]-temp.2[[2]]), na.rm = TRUE))
    cp.c_c  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[5]]-temp.1[[2]])*(temp.2[[5]]-temp.2[[2]]), na.rm = TRUE))
    cp.l_r  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[3]]-temp.1[[2]])*(temp.2[[4]]-temp.2[[2]]), na.rm = TRUE))
    cp.l_c  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[3]]-temp.1[[2]])*(temp.2[[5]]-temp.2[[2]]), na.rm = TRUE))
    cp.r_l  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[4]]-temp.1[[2]])*(temp.2[[3]]-temp.2[[2]]), na.rm = TRUE))
    cp.r_c  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[4]]-temp.1[[2]])*(temp.2[[5]]-temp.2[[2]]), na.rm = TRUE))
    cp.c_l  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[5]]-temp.1[[2]])*(temp.2[[3]]-temp.2[[2]]), na.rm = TRUE))
    cp.c_r  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[5]]-temp.1[[2]])*(temp.2[[4]]-temp.2[[2]]), na.rm = TRUE))

    cp.lr_lr  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[6]]-temp.1[[2]])*(temp.2[[6]]-temp.2[[2]]), na.rm = TRUE))
    cp.lc_lc  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[7]]-temp.1[[2]])*(temp.2[[7]]-temp.2[[2]]), na.rm = TRUE))
    cp.rc_rc  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[8]]-temp.1[[2]])*(temp.2[[8]]-temp.2[[2]]), na.rm = TRUE))
    cp.lr_rl  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[6]]-temp.1[[2]])*(temp.2[[9]]-temp.2[[2]]), na.rm = TRUE))
    cp.lc_cl  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[7]]-temp.1[[2]])*(temp.2[[10]]-temp.2[[2]]), na.rm = TRUE))
    cp.rc_cr  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[8]]-temp.1[[2]])*(temp.2[[11]]-temp.2[[2]]), na.rm = TRUE))
    cp.lr_lc  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[6]]-temp.1[[2]])*(temp.2[[7]]-temp.2[[2]]), na.rm = TRUE))
    cp.lr_cr  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[6]]-temp.1[[2]])*(temp.2[[11]]-temp.2[[2]]), na.rm = TRUE))
    cp.lc_lr  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[7]]-temp.1[[2]])*(temp.2[[6]]-temp.2[[2]]), na.rm = TRUE))
    cp.lc_rc  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[7]]-temp.1[[2]])*(temp.2[[8]]-temp.2[[2]]), na.rm = TRUE))
    cp.rc_rl  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[8]]-temp.1[[2]])*(temp.2[[9]]-temp.2[[2]]), na.rm = TRUE))
    cp.rc_lc  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[8]]-temp.1[[2]])*(temp.2[[7]]-temp.2[[2]]), na.rm = TRUE))
    cp.lr_cl  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[6]]-temp.1[[2]])*(temp.2[[10]]-temp.2[[2]]), na.rm = TRUE))
    cp.lr_rc  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[6]]-temp.1[[2]])*(temp.2[[8]]-temp.2[[2]]), na.rm = TRUE))
    cp.lc_rl  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[7]]-temp.1[[2]])*(temp.2[[9]]-temp.2[[2]]), na.rm = TRUE))
    cp.lc_cr  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[7]]-temp.1[[2]])*(temp.2[[11]]-temp.2[[2]]), na.rm = TRUE))
    cp.rc_lr  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[8]]-temp.1[[2]])*(temp.2[[6]]-temp.2[[2]]), na.rm = TRUE))
    cp.rc_cl  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[8]]-temp.1[[2]])*(temp.2[[10]]-temp.2[[2]]), na.rm = TRUE))

    cp.lrc_lrc  <- (sum((temp.1[[12]]-temp.1[[2]])*(temp.2[[12]]-temp.2[[2]]), na.rm = TRUE))
    cp.lrc_rlc  <- (sum((temp.1[[12]]-temp.1[[2]])*(temp.2[[14]]-temp.2[[2]]), na.rm = TRUE))
    cp.lrc_crl  <- (sum((temp.1[[12]]-temp.1[[2]])*(temp.2[[17]]-temp.2[[2]]), na.rm = TRUE))
    cp.lrc_lcr  <- (sum((temp.1[[12]]-temp.1[[2]])*(temp.2[[13]]-temp.2[[2]]), na.rm = TRUE))
    cp.lrc_clr  <- (sum((temp.1[[12]]-temp.1[[2]])*(temp.2[[16]]-temp.2[[2]]), na.rm = TRUE))
    cp.lrc_rcl  <- (sum((temp.1[[12]]-temp.1[[2]])*(temp.2[[15]]-temp.2[[2]]), na.rm = TRUE))

    p1  <- (temp.1[[1]]-2)*(temp.1[[1]]-1)^2
    p2  <- (temp.1[[1]]-2)
    p3  <- -(temp.1[[1]]-1)*(temp.1[[1]]-2)
    p4  <- (temp.1[[1]]-1)*(temp.1[[1]]-2)
    p5  <- -(temp.1[[1]]-2)
    p6  <- (temp.1[[1]]-1)
    p7  <- 2*(temp.1[[1]]-1)
    p8  <- (temp.1[[1]]-2)*(temp.1[[1]]^2-temp.1[[1]]-1)
    p9  <- (temp.1[[1]]^2-2*temp.1[[1]]+2)
    p10 <- (temp.1[[1]]^2-temp.1[[1]]-1)
    p11 <- (temp.1[[1]]*(temp.1[[1]]-1)*(temp.1[[1]]-2)-1)

    C <- matrix(c(p1,p2,p2,p3,p3,p3,p2,p3,p2,p4,p4,2,p5,p5,2,p4,p5,p4,p5,p5,p5,p5,p5,p5,p5,p5,p5,p6,-1,-1,p6,-1,-1,
                  p2,p1,p2,p3,p2,p3,p3,p2,p3,p4,2,p4,p5,2,p5,p5,p4,p5,p5,p4,p5,p5,p5,p5,p5,p5,p5,p6,-1,p6,-1,-1,-1,
                  p2,p2,p1,p2,p3,p2,p3,p3,p3,2,p4,p4,2,p5,p5,p5,p5,p5,p4,p5,p4,p5,p5,p5,p5,p5,p5,p6,p6,-1,-1,-1,-1,
                  p3,p3,p2,p1,p3,p2,p2,p2,p3,p5,p5,p5,p4,p5,p5,p5,p5,p5,p4,p5,2,p5,p4,p4,p5,p5,2,-1,p6,-1,-1,-1,p6,
                  p3,p2,p3,p3,p1,p2,p3,p2,p2,p5,p5,p5,p5,p4,p5,p5,p4,p5,p5,2,p5,p4,p5,p5,p4,2,p5,-1,-1,p6,-1,p6,-1,
                  p3,p3,p2,p2,p2,p1,p3,p3,p2,p5,p5,p5,p4,p5,p5,p5,p5,p5,2,p5,p4,p4,p5,p5,2,p4,p5,-1,p6,-1,-1,p6,-1,
                  p2,p3,p3,p2,p3,p3,p1,p2,p2,p5,p5,p5,p5,p5,p4,p4,p5,2,p5,p5,p5,p5,p4,2,p5,p5,p4,-1,-1,-1,p6,-1,p6,
                  p3,p2,p3,p2,p2,p3,p2,p1,p3,p5,p5,p5,p5,p4,p5,p5,2,p5,p5,p4,p5,p5,2,p4,p5,p5,p4,-1,-1,p6,-1,-1,p6,
                  p2,p3,p3,p3,p2,p2,p2,p3,p1,p5,p5,p5,p5,p5,p4,2,p5,p4,p5,p5,p5,2,p5,p5,p4,p4,p5,-1,-1,-1,p6,p6,-1,
                  p1,p1,p7,p3,p3,p3,p3,p3,p3,p8,p9,p9,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p10,-1,-1,-1,-1,-1,
                  p1,p7,p1,p3,p3,p3,p3,p3,p3,p9,p8,p9,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p10,-1,-1,-1,-1,-1,
                  p7,p1,p1,p3,p3,p3,p3,p3,p3,p9,p9,p8,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p10,-1,-1,-1,-1,-1,
                  p3,p3,p7,p1,p3,p1,p3,p3,p3,p5,p5,p5,p8,p5,p5,p5,p5,p5,p9,p5,p9,p5,p5,p5,p5,p5,p5,-1,p10,-1,-1,-1,-1,
                  p3,p7,p3,p3,p1,p3,p3,p1,p3,p5,p5,p5,p5,p8,p5,p5,p9,p5,p5,p9,p5,p5,p5,p5,p5,p5,p5,-1,-1,p10,-1,-1,-1,
                  p7,p3,p3,p3,p3,p3,p1,p3,p1,p5,p5,p5,p5,p5,p8,p9,p5,p9,p5,p5,p5,p5,p5,p5,p5,p5,p5,-1,-1,-1,p10,-1,-1,
                  p1,p3,p3,p3,p3,p3,p1,p3,p7,p5,p5,p5,p5,p5,p9,p8,p5,p9,p5,p5,p5,p5,p5,p5,p5,p5,p5,-1,-1,-1,p10,-1,-1,
                  p3,p1,p3,p3,p1,p3,p3,p7,p3,p5,p5,p5,p5,p9,p5,p5,p8,p5,p5,p9,p5,p5,p5,p5,p5,p5,p5,-1,-1,p10,-1,-1,-1,
                  p1,p3,p3,p3,p3,p3,p7,p3,p1,p5,p5,p5,p5,p5,p9,p9,p5,p8,p5,p5,p5,p5,p5,p5,p5,p5,p5,-1,-1,-1,p10,-1,-1,
                  p3,p3,p1,p1,p3,p7,p3,p3,p3,p5,p5,p5,p9,p5,p5,p5,p5,p5,p8,p5,p9,p5,p5,p5,p5,p5,p5,-1,p10,-1,-1,-1,-1,
                  p3,p1,p3,p3,p7,p3,p3,p1,p3,p5,p5,p5,p5,p9,p5,p5,p9,p5,p5,p8,p5,p5,p5,p5,p5,p5,p5,-1,-1,p10,-1,-1,-1,
                  p3,p3,p1,p7,p3,p1,p3,p3,p3,p5,p5,p5,p9,p5,p5,p5,p5,p5,p9,p5,p8,p5,p5,p5,p5,p5,p5,-1,p10,-1,-1,-1,-1,
                  p3,p3,p3,p3,p1,p1,p3,p3,p7,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p8,p5,p5,p9,p9,p5,-1,-1,-1,-1,p10,-1,
                  p3,p3,p3,p1,p3,p3,p1,p7,p3,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p8,p9,p5,p5,p9,-1,-1,-1,-1,-1,p10,
                  p3,p3,p3,p1,p3,p3,p7,p1,p3,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p9,p8,p5,p5,p9,-1,-1,-1,-1,-1,p10,
                  p3,p3,p3,p3,p1,p7,p3,p3,p1,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p9,p5,p5,p8,p9,p5,-1,-1,-1,-1,p10,-1,
                  p3,p3,p3,p3,p7,p1,p3,p3,p1,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p9,p5,p5,p9,p8,p5,-1,-1,-1,-1,p10,-1,
                  p3,p3,p3,p7,p3,p3,p1,p1,p3,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p9,p9,p5,p5,p8,-1,-1,-1,-1,-1,p10,
                  p1,p1,p1,p3,p3,p3,p3,p3,p3,p8,p8,p8,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p11,-1,-1,-1,-1,-1,
                  p3,p3,p1,p1,p3,p1,p3,p3,p3,p5,p5,p5,p8,p5,p5,p5,p5,p5,p8,p5,p8,p5,p5,p5,p5,p5,p5,-1,p11,-1,-1,-1,-1,
                  p3,p1,p3,p3,p1,p3,p3,p1,p3,p5,p5,p5,p5,p8,p5,p5,p8,p5,p5,p8,p5,p5,p5,p5,p5,p5,p5,-1,-1,p11,-1,-1,-1,
                  p1,p3,p3,p3,p3,p3,p1,p3,p1,p5,p5,p5,p5,p5,p8,p8,p5,p8,p5,p5,p5,p5,p5,p5,p5,p5,p5,-1,-1,-1,p11,-1,-1,
                  p3,p3,p3,p3,p1,p1,p3,p3,p1,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p8,p5,p5,p8,p8,p5,-1,-1,-1,-1,p11,-1,
                  p3,p3,p3,p1,p3,p3,p1,p1,p3,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p5,p8,p8,p5,p5,p8,-1,-1,-1,-1,-1,p11),
                nrow=33, ncol=33, byrow=TRUE)

    if (temp.1[[1]]>5){
      S     <- c(cp.l_l, cp.r_r, cp.c_c, cp.l_r, cp.l_c, cp.r_l, cp.r_c, cp.c_l, cp.c_r,
                 cp.lr_lr, cp.lc_lc, cp.rc_rc, cp.lr_rl, cp.lc_cl, cp.rc_cr,
                 cp.lr_lc, cp.lr_cr, cp.lc_lr, cp.lc_rc, cp.rc_rl, cp.rc_lc,
                 cp.lr_cl, cp.lr_rc, cp.lc_rl, cp.lc_cr, cp.rc_lr, cp.rc_cl,
                 cp.lrc_lrc, cp.lrc_rlc, cp.lrc_crl, cp.lrc_lcr, cp.lrc_clr, cp.lrc_rcl)

      IC    <- solve(C)
      cor.1 <- as.vector(IC%*%S)
    }else{
      S     <- c(cp.l_l, cp.r_r, cp.c_c, cp.l_r, cp.l_c, cp.r_l, cp.r_c, cp.c_l, cp.c_r,
                 cp.lr_lr, cp.lc_lc, cp.rc_rc, cp.lr_rl, cp.lc_cl, cp.rc_cr,
                 cp.lr_lc, cp.lr_cr, cp.lc_lr, cp.lc_rc, cp.rc_rl, cp.rc_lc,
                 cp.lr_cl, cp.lr_rc, cp.lc_rl, cp.lc_cr, cp.rc_lr, cp.rc_cl,
                 cp.lrc_lrc, cp.lrc_rlc, cp.lrc_crl, cp.lrc_lcr)
      C     <- C[-32:-33,-32:-33]
      IC    <- solve(C)
      cor.1 <- as.vector(IC%*%S)
      cor.1 <- c(cor.1, NA, NA)
    }
    cor.list[[i]] <- list(temp.1[[1]], cor.1)
  }
  sf  <- trmGlobal$suffixes

  temp.2  <- c("n", paste(var.1,sf[1]," with ",var.2,sf[1]," covariance", sep = ""),
               paste(var.1,sf[2]," with ",var.2,sf[2]," covariance", sep = ""),
               paste(var.1,sf[3]," with ",var.2,sf[3]," covariance", sep = ""),
               paste(var.1,sf[1]," with ",var.2,sf[2]," covariance", sep = ""),
               paste(var.1,sf[1]," with ",var.2,sf[3]," covariance", sep = ""),
               paste(var.1,sf[2]," with ",var.2,sf[1]," covariance", sep = ""),
               paste(var.1,sf[2]," with ",var.2,sf[3]," covariance", sep = ""),
               paste(var.1,sf[3]," with ",var.2,sf[1]," covariance", sep = ""),
               paste(var.1,sf[3]," with ",var.2,sf[2]," covariance", sep = ""),

               paste(var.1,sf[4]," with ",var.2,sf[4]," covariance", sep = ""),
               paste(var.1,sf[5]," with ",var.2,sf[5]," covariance", sep = ""),
               paste(var.1,sf[6]," with ",var.2,sf[6]," covariance", sep = ""),

               paste(var.1,sf[4]," with ",var.2,sf[7]," covariance", sep = ""),
               paste(var.1,sf[5]," with ",var.2,sf[8]," covariance", sep = ""),
               paste(var.1,sf[6]," with ",var.2,sf[9]," covariance", sep = ""),

               paste(var.1,sf[4]," with ",var.2,sf[5]," covariance", sep = ""),
               paste(var.1,sf[4]," with ",var.2,sf[9]," covariance", sep = ""),
               paste(var.1,sf[5]," with ",var.2,sf[4]," covariance", sep = ""),

               paste(var.1,sf[5]," with ",var.2,sf[6]," covariance", sep = ""),
               paste(var.1,sf[6]," with ",var.2,sf[7]," covariance", sep = ""),
               paste(var.1,sf[6]," with ",var.2,sf[5]," covariance", sep = ""),

               paste(var.1,sf[4]," with ",var.2,sf[8]," covariance", sep = ""),
               paste(var.1,sf[4]," with ",var.2,sf[6]," covariance", sep = ""),
               paste(var.1,sf[5]," with ",var.2,sf[7]," covariance", sep = ""),

               paste(var.1,sf[5]," with ",var.2,sf[9]," covariance", sep = ""),
               paste(var.1,sf[6]," with ",var.2,sf[4]," covariance", sep = ""),
               paste(var.1,sf[6]," with ",var.2,sf[8]," covariance", sep = ""),

               paste(var.1,sf[10]," with ",var.2,sf[10]," covariance", sep = ""),
               paste(var.1,sf[10]," with ",var.2,sf[13]," covariance", sep = ""),
               paste(var.1,sf[10]," with ",var.2,sf[15]," covariance", sep = ""),
               paste(var.1,sf[10]," with ",var.2,sf[11]," covariance", sep = ""),
               paste(var.1,sf[10]," with ",var.2,sf[14]," covariance", sep = ""),
               paste(var.1,sf[10]," with ",var.2,sf[12]," covariance", sep = ""))

  res  <- setNames(data.frame(matrix(unlist(cor.list), nrow=length(var1.grr.means), byrow=TRUE)),
                   temp.2)
  res = cbind(group = names(var1.grr.means), res)
  #print(res)
  return(res)
}

trm.univariate.cor <- function(var.grr.means, var.1){

  cor.list  <- list()
  ### Calculate correlations for each group
  for (i in names(var.grr.means)){
    temp.1  <- var.grr.means[[i]]

    ss.l    <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[3]]-temp.1[[2]])^2, na.rm = TRUE))
    ss.r    <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[4]]-temp.1[[2]])^2, na.rm = TRUE))
    ss.c    <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[5]]-temp.1[[2]])^2, na.rm = TRUE))

    ss.lr   <- ((temp.1[[1]]-2)*
                  sum((temp.1[[6]]-temp.1[[2]])^2, na.rm = TRUE))
    ss.lc   <- ((temp.1[[1]]-2)*
                  sum((temp.1[[7]]-temp.1[[2]])^2, na.rm = TRUE))
    ss.rc   <- ((temp.1[[1]]-2)*
                  sum((temp.1[[8]]-temp.1[[2]])^2, na.rm = TRUE))

    ss.lrc  <- (sum((temp.1[[12]]-temp.1[[2]])^2, na.rm = TRUE))

    cp.l_r  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[3]]-temp.1[[2]])*(temp.1[[4]]-temp.1[[2]]), na.rm = TRUE))
    cp.l_c  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[3]]-temp.1[[2]])*(temp.1[[5]]-temp.1[[2]]), na.rm = TRUE))
    cp.r_c  <- ((temp.1[[1]]-1)*(temp.1[[1]]-2)*
                  sum((temp.1[[4]]-temp.1[[2]])*(temp.1[[5]]-temp.1[[2]]), na.rm = TRUE))

    cp.lr_rl  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[6]]-temp.1[[2]])*(temp.1[[9]]-temp.1[[2]]), na.rm = TRUE))
    cp.lc_cl  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[7]]-temp.1[[2]])*(temp.1[[10]]-temp.1[[2]]), na.rm = TRUE))
    cp.rc_cr  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[8]]-temp.1[[2]])*(temp.1[[11]]-temp.1[[2]]), na.rm = TRUE))
    cp.lr_lc  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[6]]-temp.1[[2]])*(temp.1[[7]]-temp.1[[2]]), na.rm = TRUE))
    cp.lr_cr  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[6]]-temp.1[[2]])*(temp.1[[11]]-temp.1[[2]]), na.rm = TRUE))
    cp.lc_rc  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[7]]-temp.1[[2]])*(temp.1[[8]]-temp.1[[2]]), na.rm = TRUE))
    cp.lr_cl  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[6]]-temp.1[[2]])*(temp.1[[10]]-temp.1[[2]]), na.rm = TRUE))
    cp.lr_rc  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[6]]-temp.1[[2]])*(temp.1[[8]]-temp.1[[2]]), na.rm = TRUE))
    cp.lc_cr  <- ((temp.1[[1]]-2)*
                    sum((temp.1[[7]]-temp.1[[2]])*(temp.1[[11]]-temp.1[[2]]), na.rm = TRUE))

    cp.lrc_rlc  <- (sum((temp.1[[12]]-temp.1[[2]])*(temp.1[[14]]-temp.1[[2]]), na.rm = TRUE))
    cp.lrc_crl  <- (sum((temp.1[[12]]-temp.1[[2]])*(temp.1[[17]]-temp.1[[2]]), na.rm = TRUE))
    cp.lrc_lcr  <- (sum((temp.1[[12]]-temp.1[[2]])*(temp.1[[13]]-temp.1[[2]]), na.rm = TRUE))
    cp.lrc_rcl  <- (sum((temp.1[[12]]-temp.1[[2]])*((temp.1[[15]]+temp.1[[16]])/2-temp.1[[2]]), na.rm = TRUE))

    p1  <- (temp.1[[1]]-2)*(temp.1[[1]]-1)^2
    p2  <- (temp.1[[1]]-2)
    p3  <- (temp.1[[1]]-1)*(temp.1[[1]]-2)
    p4  <- (temp.1[[1]]-1)
    p5  <- -2*(temp.1[[1]]-1)*(temp.1[[1]]-2)
    p6  <- -1*(temp.1[[1]]-2)
    p7  <- 2*(temp.1[[1]]-2)
    p8  <- -2*(temp.1[[1]]-2)
    p9  <- 2*(temp.1[[1]]-1)*(temp.1[[1]]-2)
    p10 <- 2*(temp.1[[1]]-1)
    p11 <- (temp.1[[1]]-2)*(temp.1[[1]]^2-temp.1[[1]]-1)
    p12 <- (temp.1[[1]]^2-2*temp.1[[1]]+2)
    p13 <- (temp.1[[1]]^2-temp.1[[1]]-1)
    p14 <- temp.1[[1]]*(temp.1[[1]]-1)*(temp.1[[1]]-2)-1
    p15 <- -1*(temp.1[[1]]-1)*(temp.1[[1]]-2)
    p16 <- (temp.1[[1]]-2)*(temp.1[[1]]^2-2*temp.1[[1]]+2)
    p17 <- -1*(temp.1[[1]]-2)^2
    p18 <- -1*(temp.1[[1]]-4)
    p19 <- (temp.1[[1]]-2)^2
    p20 <- (temp.1[[1]]-1)*(temp.1[[1]]-2)+2
    p21 <- 2*(temp.1[[1]]-2)*(temp.1[[1]]-1)^2
    p22 <- (temp.1[[1]]+1)*(temp.1[[1]]-1)*(temp.1[[1]]-2)+2
    p23 <- 2*(temp.1[[1]]^2-2*temp.1[[1]]+2)
    p24 <- 2*(temp.1[[1]]-1)+(temp.1[[1]]-2)*(temp.1[[1]]-1)^2
    p25 <- (temp.1[[1]]-1)*(temp.1[[1]]-2)^2
    p26 <- (temp.1[[1]]-2)*(temp.1[[1]]^2-temp.1[[1]]-2)
    p27 <- -1*(temp.1[[1]]-1)*(temp.1[[1]]-4)
    p28 <- temp.1[[1]]+(temp.1[[1]]-2)^2
    p29 <- temp.1[[1]]*(temp.1[[1]]-1)-2
    p30 <- 2*(temp.1[[1]]-2)*(temp.1[[1]]^2-temp.1[[1]]-1)
    p31 <- temp.1[[1]]*(temp.1[[1]]-1)*(temp.1[[1]]-2)-2

    C <- matrix(c(p1,p2,p2,p3,p3,2,p4,p5,p5,p7,p6,p6,2,p9,p8,p8,p8,p8,p8,-1,-1,p4,-2,
                  p2,p1,p2,p3,2,p3,p4,p5,p7,p5,p6,2,p6,p8,p9,p8,p8,p8,p8,-1,p4,-1,-2,
                  p2,p2,p1,2,p3,p3,p4,p7,p5,p5,2,p6,p6,p8,p8,p9,p8,p8,p8,p4,-1,-1,-2,
                  p1,p1,p10,p11,p12,p12,p13,p5,p5,p5,p6,p6,p6,p8,p8,p8,p8,p8,p8,-1,-1,-1,-2,
                  p1,p10,p1,p12,p11,p12,p13,p5,p5,p5,p6,p6,p6,p8,p8,p8,p8,p8,p8,-1,-1,-1,-2,
                  p10,p1,p1,p12,p12,p11,p13,p5,p5,p5,p6,p6,p6,p8,p8,p8,p8,p8,p8,-1,-1,-1,-2,
                  p1,p1,p1,p11,p11,p11,p14,p5,p5,p5,p6,p6,p6,p8,p8,p8,p8,p8,p8,-1,-1,-1,-2,
                  p15,p15,p2,p6,p6,p6,-1,p16,p17,p17,p3,p6,p6,p8,p8,p20,p19,p19,p18,p4,-1,-1,p2,
                  p15,p2,p15,p6,p6,p6,-1,p17,p16,p17,p6,p3,p6,p8,p20,p8,p19,p18,p19,-1,p4,-1,p2,
                  p2,p15,p15,p6,p6,p6,-1,p17,p17,p16,p6,p6,p3,p20,p8,p8,p18,p19,p19,-1,-1,p4,p2,
                  p15,p15,p10,p6,p6,p6,-1,p21,p5,p5,p11,p6,p6,p8,p8,p23,p8,p8,p8,p13,-1,-1,-2,
                  p15,p10,p15,p6,p6,p6,-1,p5,p21,p5,p6,p11,p6,p8,p23,p8,p8,p8,p8,-1,p13,-1,-2,
                  p10,p15,p15,p6,p6,p6,-1,p5,p5,p21,p6,p6,p11,p23,p8,p8,p8,p8,p8,-1,-1,p13,-2,
                  p1,p15,p15,p6,p6,p6,-1,p5,p5,p24,p6,p6,p12,p22,p8,p8,p8,p8,p8,-1,-1,p13,-2,
                  p15,p1,p15,p6,p6,p6,-1,p5,p24,p5,p6,p12,p6,p8,p22,p8,p8,p8,p8,-1,p13,-1,-2,
                  p15,p15,p1,p6,p6,p6,-1,p24,p5,p5,p12,p6,p6,p8,p8,p22,p8,p8,p8,p13,-1,-1,-2,
                  p15,p15,p15,p6,p6,p6,-1,p25,p25,p27,p6,p6,p6,p8,p8,p8,p26,p28,p28,-1,-1,-1,p29,
                  p15,p15,p15,p6,p6,p6,-1,p25,p27,p25,p6,p6,p6,p8,p8,p8,p28,p26,p28,-1,-1,-1,p29,
                  p15,p15,p15,p6,p6,p6,-1,p27,p25,p25,p6,p6,p6,p8,p8,p8,p28,p28,p26,-1,-1,-1,p29,
                  p15,p15,p1,p6,p6,p6,-1,p21,p5,p5,p11,p6,p6,p8,p8,p30,p8,p8,p8,p14,-1,-1,-2,
                  p15,p1,p15,p6,p6,p6,-1,p5,p21,p5,p6,p11,p6,p8,p30,p8,p8,p8,p8,-1,p14,-1,-2,
                  p1,p15,p15,p6,p6,p6,-1,p5,p5,p21,p6,p6,p11,p30,p8,p8,p8,p8,p8,-1,-1,p14,-2,
                  p15,p15,p15,p6,p6,p6,-1,p25,p25,p25,p6,p6,p6,p8,p8,p8,p26,p26,p26,-1,-1,-1,p31),
                nrow=23, ncol=23, byrow=TRUE)

    if (temp.1[[1]]>5){
      S     <- c(ss.l, ss.r, ss.c, ss.lr, ss.lc, ss.rc, ss.lrc,
                 cp.l_r, cp.l_c, cp.r_c, cp.lr_rl, cp.lc_cl, cp.rc_cr,
                 cp.lr_lc, cp.lr_cr, cp.lc_rc, cp.lr_cl, cp.lr_rc, cp.lc_cr,
                 cp.lrc_rlc, cp.lrc_crl, cp.lrc_lcr, cp.lrc_rcl)
      IC    <- solve(C)
      cor.1 <- as.vector(IC%*%S)
    }else{
      S     <- c(ss.l, ss.r, ss.c, ss.lr, ss.lc, ss.rc, ss.lrc,
                 cp.l_r, cp.l_c, cp.r_c, cp.lr_rl, cp.lc_cl, cp.rc_cr,
                 cp.lr_lc, cp.lr_cr, cp.lc_rc, cp.lr_cl, cp.lr_rc, cp.lc_cr,
                 cp.lrc_rlc, cp.lrc_crl, cp.lrc_lcr)
      C     <- C[-23,-23]
      IC    <- solve(C)
      cor.1 <- as.vector(IC%*%S)
      cor.1 <- c(cor.1, NA)
    }
    ####
    cor.list[[i]] <- list(temp.1[[1]], temp.1[[2]], cor = cor.1)
  }
  #####
  sf  <- trmGlobal$suffixes

  temp.2  <- c("n", "group.mean",
               paste(var.1,sf[1]," variance", sep = ""),
               paste(var.1,sf[2]," variance", sep = ""),
               paste(var.1,sf[3]," variance", sep = ""),

               paste(var.1,sf[4]," variance", sep = ""),
               paste(var.1,sf[5]," variance", sep = ""),
               paste(var.1,sf[6]," variance", sep = ""),

               paste(var.1,sf[10]," variance", sep = ""),

               paste(var.1,sf[1]," with ",var.1,sf[2]," covariance", sep = ""),
               paste(var.1,sf[1]," with ",var.1,sf[3]," covariance", sep = ""),
               paste(var.1,sf[2]," with ",var.1,sf[3]," covariance", sep = ""),

               paste(var.1,sf[4]," with ",var.1,sf[7]," covariance", sep = ""),
               paste(var.1,sf[5]," with ",var.1,sf[8]," covariance", sep = ""),
               paste(var.1,sf[6]," with ",var.1,sf[9]," covariance", sep = ""),

               paste(var.1,sf[4]," with ",var.1,sf[5]," covariance", sep = ""),
               paste(var.1,sf[4]," with ",var.1,sf[9]," covariance", sep = ""),
               paste(var.1,sf[5]," with ",var.1,sf[6]," covariance", sep = ""),

               paste(var.1,sf[4]," with ",var.1,sf[8]," covariance", sep = ""),
               paste(var.1,sf[4]," with ",var.1,sf[6]," covariance", sep = ""),
               paste(var.1,sf[5]," with ",var.1,sf[9]," covariance", sep = ""),

               paste(var.1,sf[10]," with ",var.1,sf[13]," covariance", sep = ""),
               paste(var.1,sf[10]," with ",var.1,sf[15]," covariance", sep = ""),
               paste(var.1,sf[10]," with ",var.1,sf[11]," covariance", sep = ""),
               paste(var.1,sf[10]," with ",var.1,sf[12]," covariance", sep = ""))


  res  <- setNames(data.frame(matrix(unlist(cor.list), nrow=length(var.grr.means), byrow=TRUE)),
                   temp.2)
  res = cbind(group = names(var.grr.means), res)
  #print(res)
  return(res)
}


trm.effects <- function(grr.means, var.1, group.id, layer.id, row.id, column.id){

  ind.effects.l  <- data.frame()
  ind.effects.r  <- data.frame()
  ind.effects.c  <- data.frame()
  dya.effects.lr <- data.frame()
  dya.effects.lc <- data.frame()
  dya.effects.rc <- data.frame()
  tri.effects.lrc <- data.frame()

  ### Calculate effects for each group
  for (i in names(grr.means)){
    temp.1      <- grr.means[[i]]
    n           <- temp.1[[1]]

    g.means     <- temp.1[[2]]

    l.means     <- temp.1[[3]]
    r.means     <- temp.1[[4]]
    c.means     <- temp.1[[5]]

    lr.means    <- temp.1[[6]]
    lc.means    <- temp.1[[7]]
    rc.means    <- temp.1[[8]]
    rl.means    <- temp.1[[9]]
    cl.means    <- temp.1[[10]]
    cr.means    <- temp.1[[11]]

    l.effects   <- ((n-1)*(n-2)/(n*(n-3))*l.means
                    +(n-1)/(n*(n-3))*r.means
                    +(n-1)/(n*(n-3))*c.means
                    -(n-1)/(n-3)*g.means
                    )             ### Judge Effects
    r.effects   <- ((n-1)*(n-2)/(n*(n-3))*r.means
                    +(n-1)/(n*(n-3))*l.means
                    +(n-1)/(n*(n-3))*c.means
                    -(n-1)/(n-3)*g.means
                    )             ### Actor Effects
    c.effects   <- ((n-1)*(n-2)/(n*(n-3))*c.means
                    +(n-1)/(n*(n-3))*l.means
                    +(n-1)/(n*(n-3))*r.means
                    -(n-1)/(n-3)*g.means
                    )             ### Partner Effects

    l.mat2d     <- array(l.effects,c(n,n))
    r.mat2d     <- array(r.effects,c(n,n))
    c.mat2d     <- array(c.effects,c(n,n))

    lr.effects  <- ((n-2)^2*(n^2-4*n+1)/(n*(n-1)*(n-3)*(n-4))*lr.means
                    +2*(n-2)/(n*(n-1)*(n-3)*(n-4))*rl.means
                    +(n-2)*(n^2-4*n+2)/(n*(n-1)*(n-3)*(n-4))*(cr.means+lc.means)
                    +(n-2)^2/(n*(n-1)*(n-3)*(n-4))*(rc.means+cl.means)
                    -(n^3-7*n^2+14*n-7)/((n-1)*(n-3)*(n-4))*(l.mat2d+t(r.mat2d))
                    -((n-1)*(n-3)*(n-4))^(-1)*(r.mat2d+t(l.mat2d))
                    -((n-1)*(n-4))^(-1)*(c.mat2d+t(c.mat2d))
                    -(n-2)/(n-4)*g.means
                    )             ### Judge-Actor Effects
    lc.effects  <- ((n-2)^2*(n^2-4*n+1)/(n*(n-1)*(n-3)*(n-4))*lc.means
                    +2*(n-2)/(n*(n-1)*(n-3)*(n-4))*cl.means
                    +(n-2)*(n^2-4*n+2)/(n*(n-1)*(n-3)*(n-4))*(rc.means+lr.means)
                    +(n-2)^2/(n*(n-1)*(n-3)*(n-4))*(cr.means+rl.means)
                    -(n^3-7*n^2+14*n-7)/((n-1)*(n-3)*(n-4))*(l.mat2d+t(c.mat2d))
                    -((n-1)*(n-3)*(n-4))^(-1)*(c.mat2d+t(l.mat2d))
                    -((n-1)*(n-4))^(-1)*(r.mat2d+t(r.mat2d))
                    -(n-2)/(n-4)*g.means
                    )             ### Judge-Partner Effects
    rc.effects  <- ((n-2)^2*(n^2-4*n+1)/(n*(n-1)*(n-3)*(n-4))*rc.means
                    +2*(n-2)/(n*(n-1)*(n-3)*(n-4))*cr.means
                    +(n-2)*(n^2-4*n+2)/(n*(n-1)*(n-3)*(n-4))*(rl.means+lc.means)
                    +(n-2)^2/(n*(n-1)*(n-3)*(n-4))*(lr.means+cl.means)
                    -(n^3-7*n^2+14*n-7)/((n-1)*(n-3)*(n-4))*(r.mat2d+t(c.mat2d))
                    -((n-1)*(n-3)*(n-4))^(-1)*(c.mat2d+t(r.mat2d))
                    -((n-1)*(n-4))^(-1)*(l.mat2d+t(l.mat2d))
                    -(n-2)/(n-4)*g.means
                    )             ### Actor-Partner Effects
    #rl.effects  <- t(lr.effects)  ### Actor-Judge Effects
    #cl.effects  <- t(lc.effects)  ### Partner-Judge Effects
    #cr.effects  <- t(rc.effects)  ### Partner-Actor Effects

    l.mat3d     <- array(l.effects,c(n,n,n))
    r.mat3d     <- array(r.effects,c(n,n,n))
    c.mat3d     <- array(c.effects,c(n,n,n))
    lr.mat3d    <- array(lr.effects,c(n,n,n))
    lc.mat3d    <- array(lc.effects,c(n,n,n))
    rc.mat3d    <- array(rc.effects,c(n,n,n))

    lrc.effects <- (temp.1[[12]]
                    -aperm(l.mat3d, c(3,2,1))-aperm(r.mat3d, c(1,2,3))-aperm(c.mat3d, c(2,1,3))
                    -aperm(lr.mat3d, c(2,3,1))-aperm(lc.mat3d, c(3,2,1))-aperm(rc.mat3d, c(1,2,3))
                    -g.means)                   ### Judge-Actor-Partner Effects
    #lcr.effects <- aperm(lrc.effects, c(2,1,3)) ### Judge-Partner-Actor Effects
    #rlc.effects <- aperm(lrc.effects, c(3,2,1)) ### Actor-Judge-Partner Effects
    #rcl.effects <- aperm(lrc.effects, c(2,3,1)) ### Actor-Partner-Judge Effects
    #clr.effects <- aperm(lrc.effects, c(3,1,2)) ### Partner-Judge-Actor Effects
    #crl.effects <- aperm(lrc.effects, c(1,3,2)) ### Partner-Actor-Judge Effects

    #a=as.data.frame.table(l.effects)

    temp.2  <- cbind(i, id = names(temp.1[[3]]))
    l.e.long  <- cbind(temp.2, as.data.frame(l.effects))
    r.e.long  <- cbind(temp.2, as.data.frame(r.effects))
    c.e.long  <- cbind(temp.2, as.data.frame(c.effects))

    ind.effects.l  <- rbind(ind.effects.l, l.e.long)
    ind.effects.r  <- rbind(ind.effects.r, r.e.long)
    ind.effects.c  <- rbind(ind.effects.c, c.e.long)

    lr.e.long   <- cbind(i, as.data.frame.table(lr.effects))
    lc.e.long   <- cbind(i, as.data.frame.table(lc.effects))
    rc.e.long   <- cbind(i, as.data.frame.table(rc.effects))
    lrc.e.long  <- cbind(i, as.data.frame.table(lrc.effects))

    dya.effects.lr  <- rbind(dya.effects.lr, lr.e.long)
    dya.effects.lc  <- rbind(dya.effects.lc, lc.e.long)
    dya.effects.rc  <- rbind(dya.effects.rc, rc.e.long)
    tri.effects.lrc <- rbind(tri.effects.lrc, lrc.e.long)
  }

  sf  <- trmGlobal$suffixes

  ind.effects.l <- ind.effects.l[!is.na(ind.effects.l[,3]),]
  rownames(ind.effects.l) <- NULL
  colnames(ind.effects.l) <- c(group.id, layer.id, paste(var.1,sf[1], sep = ""))
  #print(ind.effects.l)

  ind.effects.r <- ind.effects.r[!is.na(ind.effects.r[,3]),]
  rownames(ind.effects.r) <- NULL
  colnames(ind.effects.r) <- c(group.id, row.id, paste(var.1,sf[2], sep = ""))
  #print(ind.effects.r)

  ind.effects.c <- ind.effects.c[!is.na(ind.effects.c[,3]),]
  rownames(ind.effects.c) <- NULL
  colnames(ind.effects.c) <- c(group.id, column.id, paste(var.1,sf[3], sep = ""))
  #print(ind.effects.c)

  dya.effects.lr <- dya.effects.lr[!is.na(dya.effects.lr[,4]),]
  rownames(dya.effects.lr) <- NULL
  colnames(dya.effects.lr) <- c(group.id, row.id, layer.id, paste(var.1,sf[4], sep = ""))
  #print(dya.effects.lr)

  dya.effects.lc <- dya.effects.lc[!is.na(dya.effects.lc[,4]),]
  rownames(dya.effects.lc) <- NULL
  colnames(dya.effects.lc) <- c(group.id, layer.id, column.id, paste(var.1,sf[5], sep = ""))
  #print(dya.effects.lc)
  #print("yo")

  dya.effects.rc <- dya.effects.rc[!is.na(dya.effects.rc[,4]),]
  rownames(dya.effects.rc) <- NULL
  colnames(dya.effects.rc) <- c(group.id, row.id, column.id, paste(var.1,sf[6], sep = ""))
  #print(dya.effects.rc)


  tri.effects.lrc <- tri.effects.lrc[!is.na(tri.effects.lrc[,5]),]
  rownames(tri.effects.lrc) <- NULL
  colnames(tri.effects.lrc) <- c(group.id, row.id, column.id, layer.id, paste(var.1,sf[10], sep = ""))
  #print(tri.effects.lrc)

  res <- list(ind.effects.l,
              ind.effects.r,
              ind.effects.c,
              dya.effects.lr[,c(1,3,2,4)],
              dya.effects.lc,
              dya.effects.rc,
              tri.effects.lrc[,c(1,4,2,3,5)])
  return(res)
}

trm.grr.means <- function(GRR){
  res <- list()
  ### Calculate effects for each group
  for (i in names(GRR)){
    temp.1  <- GRR[[i]]
    temp.2  <- list(n = nrow(temp.1),
                    g.means = mean(temp.1, na.rm = TRUE),
                    l.means = apply(temp.1, c(3), mean, na.rm = TRUE),
                    r.means = apply(temp.1, c(1), mean, na.rm = TRUE),
                    c.means = apply(temp.1, c(2), mean, na.rm = TRUE),
                    lr.means = apply(temp.1, c(3,1), mean, na.rm = TRUE),
                    lc.means = apply(temp.1, c(3,2), mean, na.rm = TRUE),
                    rc.means = apply(temp.1, c(1,2), mean, na.rm = TRUE),
                    rl.means = t(apply(temp.1, c(3,1), mean, na.rm = TRUE)),
                    cl.means = t(apply(temp.1, c(3,2), mean, na.rm = TRUE)),
                    cr.means = t(apply(temp.1, c(1,2), mean, na.rm = TRUE)),
                    var.lrc = temp.1,
                    var.lcr = aperm(temp.1, c(2,1,3)),
                    var.rlc = aperm(temp.1, c(3,2,1)),
                    var.rcl = aperm(temp.1, c(2,3,1)),
                    var.clr = aperm(temp.1, c(3,1,2)),
                    var.crl = aperm(temp.1, c(1,3,2))
    )
    res[[i]]  <- temp.2
  }
  return(res)
}

trm.create.grr <- function(trm.check, data){

  data <- data.frame(data)  # make sure that the format of the dataset is in dataframe

  if (!is.na(match(trm.check[[1]], colnames(data)))){ # Check the presence of layer.id
    if (!is.na(match(trm.check[[2]], colnames(data)))){ # Check the presence of row.id
      if (!is.na(match(trm.check[[3]], colnames(data)))){ # Check the presence of column.id
        if (trm.check[[8]] == TRUE){  # Check if is.group.id equals TRUE
          if (!is.na(match(trm.check[[4]], colnames(data)))){ # Check the presence of group.id
            if (!is.na(match(TRUE,duplicated(data[,c(trm.check[[1]],
                                                     trm.check[[2]],
                                                     trm.check[[3]],
                                                     trm.check[[4]])])))){ # check if there are duplicated record
              stop("The dataset contains duplicated record, please check.")
            }
            else{
              if (!is.na(match(trm.check[[5]], colnames(data)))){  # Check the presence of var.1
                if (!is.numeric(data[,trm.check[[5]]])){  # Check if var.1 is numeric
                  stop(trm.check[[5]], " contains non-nurmeric values.")
                }
                else{
                  if (trm.check[[7]] == FALSE){  # Check if is.univar equals FALSE
                    if (!is.na(match(trm.check[[6]], colnames(data)))){  # Check the presence of var.2
                      if (!is.numeric(data[,trm.check[[6]]])){  # Check if var.2 is numeric
                        stop(trm.check[[6]], " contains non-nurmeric values.")
                      }
                    }
                    else{ stop("Cannot find ", trm.check[[6]], " in the dataset.")}
                  }
                }
              }
              else{ stop("Cannot find ", trm.check[[5]], " in the dataset.")}
            }
          }
          else{ stop("Cannot find ", trm.check[[4]], " in the dataset.")}
        }
        else{
          if (!is.na(match(TRUE,duplicated(data[,c(trm.check[[1]],
                                                   trm.check[[2]],
                                                   trm.check[[3]])])))){ # check if there are duplicated record
            stop("The dataset contains duplicated record, please check.")
          }
          else{
            if (!is.na(match(trm.check[[5]], colnames(data)))){  # Check the presence of var.1
              if (!is.numeric(data[,trm.check[[5]]])){  # Check if var.1 is numeric
                stop(trm.check[[5]], " contains non-nurmeric values.")
              }
              else{
                if (trm.check[[7]] == FALSE){  # Check if is.univar equals FALSE
                  if (!is.na(match(trm.check[[6]], colnames(data)))){  # Check the presence of var.2
                    if (!is.numeric(data[,trm.check[[6]]])){  # Check if var.2 is numeric
                      stop(trm.check[[6]], " contains non-nurmeric values.")
                    }
                  }
                  else{ stop("Cannot find ", trm.check[[6]], " in the dataset.")}
                }
              }
            }
            else{ stop("Cannot find ", trm.check[[5]], " in the dataset.")}
          }
        }
      }
      else{ stop("Cannot find ", trm.check[[3]], " in the dataset.")}
    }
    else{ stop("Cannot find ", trm.check[[2]], " in the dataset.")}
  }
  else{ stop("Cannot find ", trm.check[[1]], " in the dataset.")}

  #### Sort the order
  data <- data[order(data[,trm.check[[3]]]),]
  data <- data[order(data[,trm.check[[2]]]),]
  data <- data[order(data[,trm.check[[1]]]),]
  if (trm.check[[8]] == TRUE){
    data <- data[order(data[,trm.check[[4]]]),]
  }

  var1.mat <- list()
  var2.mat <- list()

  ### Remove missing
  if (trm.check[[7]] == TRUE){
    data.var1 <- data[!is.na(data[,trm.check[[5]]]),]
  }else{
    data.var1 <- data[!is.na(data[,trm.check[[5]]]) | !is.na(data[,trm.check[[6]]]),]
    data.var2 <- data[!is.na(data[,trm.check[[5]]]) | !is.na(data[,trm.check[[6]]]),]
  }

  if (trm.check[[8]] == TRUE){
    group.vector <- data.var1[!duplicated(data.var1[,trm.check[[4]]]),trm.check[[4]]]
    for (i in group.vector){
      member    <- unique(c(data.var1[data.var1[,trm.check[[4]]] == i,trm.check[[1]]],
                         data.var1[data.var1[,trm.check[[4]]] == i,trm.check[[2]]],
                         data.var1[data.var1[,trm.check[[4]]] == i,trm.check[[3]]]))
      # remove groups with less than 5 members
      # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      # implement later

      temp.1    <- array(NA,c(length(member),length(member),length(member)),
                         dimnames = list(member, member, member))
      temp.2    <- data.var1[data.var1[,trm.check[[4]]] == i &
                               !(data.var1[,trm.check[[1]]] == data.var1[,trm.check[[2]]] |
                                   data.var1[,trm.check[[1]]] == data.var1[,trm.check[[3]]] |
                                   data.var1[,trm.check[[2]]] == data.var1[,trm.check[[3]]]), ]
      for (j in 1:nrow(temp.2)){
        temp.1[temp.2[j,trm.check[[2]]],
               temp.2[j,trm.check[[3]]],
               temp.2[j,trm.check[[1]]]
               ]  <- temp.2[j,trm.check[[5]]]
      }
      var1.mat[[i]]  <- temp.1
    }
  }else{
    group.vector <- c(1)
    member    <- unique(c(data.var1[ ,trm.check[[1]]],
                          data.var1[ ,trm.check[[2]]],
                          data.var1[ ,trm.check[[3]]]))

    # remove groups with less than 5 members
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # implement later

    temp.1    <- array(NA,c(length(member),length(member),length(member)),
                       dimnames = list(member, member, member))
    temp.2    <- data.var1[!(data.var1[,trm.check[[1]]] == data.var1[,trm.check[[2]]] |
                                 data.var1[,trm.check[[1]]] == data.var1[,trm.check[[3]]] |
                                 data.var1[,trm.check[[2]]] == data.var1[,trm.check[[3]]]), ]
    for (j in 1:nrow(temp.2)){
      temp.1[temp.2[j,trm.check[[2]]],
             temp.2[j,trm.check[[3]]],
             temp.2[j,trm.check[[1]]]
             ]  <- temp.2[j,trm.check[[5]]]
    }
    var1.mat[[1]]  <- temp.1
  }
  names(var1.mat) <- group.vector

  #### bivariate
  if (trm.check[[7]] == FALSE){
    if (trm.check[[8]] == TRUE){
      group.vector <- data.var2[!duplicated(data.var2[,trm.check[[4]]]),trm.check[[4]]]
      for (i in group.vector){
        member    <- unique(c(data.var2[data.var2[,trm.check[[4]]] == i,trm.check[[1]]],
                              data.var2[data.var2[,trm.check[[4]]] == i,trm.check[[2]]],
                              data.var2[data.var2[,trm.check[[4]]] == i,trm.check[[3]]]))
        # remove groups with less than 5 members
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # implement later

        temp.1    <- array(NA,c(length(member),length(member),length(member)),
                           dimnames = list(member, member, member))
        temp.2    <- data.var2[data.var2[,trm.check[[4]]] == i &
                                 !(data.var2[,trm.check[[1]]] == data.var2[,trm.check[[2]]] |
                                     data.var2[,trm.check[[1]]] == data.var2[,trm.check[[3]]] |
                                     data.var2[,trm.check[[2]]] == data.var2[,trm.check[[3]]]), ]
        for (j in 1:nrow(temp.2)){
          temp.1[temp.2[j,trm.check[[2]]],
                 temp.2[j,trm.check[[3]]],
                 temp.2[j,trm.check[[1]]]
                 ]  <- temp.2[j,trm.check[[6]]]
        }
        var2.mat[[i]]  <- temp.1
      }
    }else{
      group.vector <- c(1)
      member    <- unique(c(data.var2[ ,trm.check[[1]]],
                            data.var2[ ,trm.check[[2]]],
                            data.var2[ ,trm.check[[3]]]))

      # remove groups with less than 5 members
      # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      # implement later

      temp.1    <- array(NA,c(length(member),length(member),length(member)),
                         dimnames = list(member, member, member))
      temp.2    <- data.var2[!(data.var2[,trm.check[[1]]] == data.var2[,trm.check[[2]]] |
                                   data.var2[,trm.check[[1]]] == data.var2[,trm.check[[3]]] |
                                   data.var2[,trm.check[[2]]] == data.var2[,trm.check[[3]]]), ]
      for (j in 1:nrow(temp.2)){
        temp.1[temp.2[j,trm.check[[2]]],
               temp.2[j,trm.check[[3]]],
               temp.2[j,trm.check[[1]]]
               ]  <- temp.2[j,trm.check[[6]]]
      }
      var2.mat[[1]]  <- temp.1
    }
    names(var2.mat) <- group.vector
  }else{
    var2.mat  <- NULL
  }

  GRR  <- list(var1.mat, var2.mat)
  return(GRR)
}

trm.check <- function(formula){

  # Example:  model <- trm(var.1 ~ layer.id * row.id * column.id, data)
  # Example:  model <- trm(var.1 ~ layer.id * row.id * column.id | group.id, data)
  # Example:  model <- trm(var.1 + var.2 ~ layer.id * row.id * column.id, data)
  # Example:  model <- trm(var.1 + var.2 ~ layer.id * row.id * column.id | group.id, data)

  layer.id     <- ""    # Judge
  row.id       <- ""    # Actor
  column.id    <- ""    # Partner
  group.id     <- ""    # Group
  var.1        <- ""    #
  var.2        <- ""    #
  is.univar    <- TRUE  #
  is.group.id  <- TRUE  #

  formula <- as.formula(formula) # make sure that the input is in the formula format

  if (length(formula[[3]]) != 1){
    if (formula[[3]][[1]] == "|" && length(formula[[3]][[3]]) == 1){

      # group id is specified.
      if (length(formula[[3]][[2]]) != 1){
        if (formula[[3]][[2]][[1]] == "*" && length(formula[[3]][[2]][[2]]) != 1){
          if(formula[[3]][[2]][[2]][[1]] == "*" && length(formula[[3]][[2]][[2]][[2]]) == 1){
            if(length(formula[[2]]) == 1){

              # univariate analysis, group id is specified
              layer.id     <- formula[[3]][[2]][[2]][[2]]
              row.id       <- formula[[3]][[2]][[2]][[3]]
              column.id    <- formula[[3]][[2]][[3]]
              group.id     <- formula[[3]][[3]]
              var.1        <- formula[[2]]
            }
            else if (formula[[2]][[1]] == "+"
                      && length(formula[[2]][[2]]) == 1){

              # bivariate analysis, group id is specified
              layer.id     <- formula[[3]][[2]][[2]][[2]]
              row.id       <- formula[[3]][[2]][[2]][[3]]
              column.id    <- formula[[3]][[2]][[3]]
              group.id     <- formula[[3]][[3]]
              var.1        <- formula[[2]][[2]]
              var.2        <- formula[[2]][[3]]
              is.univar    <- FALSE
            }
            else{ stop("Wrong equation.")}
          }
          else{ stop("Wrong equation.")}
        }
        else{ stop("Wrong equation.")}
      }
      else{ stop("Wrong equation.")}
    }
    else if (formula[[3]][[1]] == "*" && length(formula[[3]][[2]]) != 1){

      # group id is not specified.
      if (formula[[3]][[2]][[1]] == "*" && length(formula[[3]][[2]][[2]]) == 1){
        if(length(formula[[2]]) == 1){

          # univariate analysis, group id is not specified
          layer.id     <- formula[[3]][[2]][[2]]
          row.id       <- formula[[3]][[2]][[3]]
          column.id    <- formula[[3]][[3]]
          group.id     <- "group.id"
          var.1        <- formula[[2]]
          is.group.id  <- FALSE
        }
        else if (formula[[2]][[1]] == "+"
                  && length(formula[[2]][[2]]) == 1){

          # bivariate analysis, group id is not specified
          layer.id     <- formula[[3]][[2]][[2]]
          row.id       <- formula[[3]][[2]][[3]]
          column.id    <- formula[[3]][[3]]
          group.id     <- "group.id"
          var.1        <- formula[[2]][[2]]
          var.2        <- formula[[2]][[3]]
          is.univar    <- FALSE
          is.group.id  <- FALSE
        }
        else{ stop("Wrong equation.")}
      }
      else{ stop("Wrong equation.")}
    }
    else{ stop("Wrong equation.")}
  }
  else{ stop("Wrong equation.")}

  res <- list(  layer.id     <- as.character(layer.id),
                row.id       <- as.character(row.id),
                column.id    <- as.character(column.id),
                group.id     <- as.character(group.id),
                var.1        <- as.character(var.1),
                var.2        <- as.character(var.2),
                is.univar    <- is.univar,
                is.group.id  <- is.group.id)

  return(res)

}
