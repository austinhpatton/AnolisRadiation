gls.ancova <- 
  function (form, Sigma, ReducedModel, FullModel, data) 
  {
    Y <- c(data[, which(colnames(data) == all.vars(form)[1])])
    n <- length(Y)
    tr <- sum(diag(Sigma))
    Sigma <- n * Sigma/tr
    invSigma <- solve(Sigma)
    df_ReducedModel <- c(length(ReducedModel[1, ]))
    df_FullModel <- c(length(FullModel[1, ]))
    C <- solve(t(ReducedModel) %*% invSigma %*% ReducedModel)
    B <- C %*% t(ReducedModel) %*% invSigma %*% Y
    SS_unexpl_ReducedModel <- (t(Y - (ReducedModel %*% B)) %*% 
                                 invSigma %*% (Y - (ReducedModel %*% B)))
    MS_unexpl_ReducedModel <- SS_unexpl_ReducedModel/c(n - df_ReducedModel)
    C <- solve(t(FullModel) %*% invSigma %*% FullModel)
    B <- C %*% t(FullModel) %*% invSigma %*% Y
    SS_unexpl_FullModel <- (t(Y - (FullModel %*% B)) %*% invSigma %*% 
                              (Y - (FullModel %*% B)))
    MS_unexpl_FullModel <- SS_unexpl_FullModel/c(n - df_FullModel)
    Fs <- ((SS_unexpl_ReducedModel - SS_unexpl_FullModel)/(c(n - 
                                                               df_ReducedModel) - c(n - df_FullModel)))/(SS_unexpl_FullModel/c(n - 
                                                                                                                                 df_FullModel))
    P <- 1 - pf(Fs, (c(n - df_ReducedModel) - c(n - df_FullModel)), 
                c(n - df_FullModel))
    Output <- as.data.frame(cbind(rbind(df_FullModel, df_ReducedModel), 
                                  rbind(round(SS_unexpl_FullModel, 4), round(SS_unexpl_ReducedModel, 
                                                                             4)), rbind(round(MS_unexpl_FullModel, 4), round(MS_unexpl_ReducedModel, 
                                                                                                                             4)), rbind(round(Fs, 4), ""), rbind(round(P, 4), 
                                                                                                                                                                 "")))
    colnames(Output) <- c("df", "Sum Sq", "Mean Sum Sq", "F value", 
                          "Pr(>F)")
    rownames(Output) <- c("FullModel", "ReducedModel")
    return(Output)
  }
