#' @export
boot632 <- function(x, y, funcs, alpha = .1, nboot = 500,  ...) {
  phat<-matrix(table(y)/length(y),ncol=1)

  theta.fit <- funcs$fitter
  theta.predict <- funcs$predictor
  err.meas <- funcs$loss

  calcse <- T
  calce00 <- T

  sein.e0<-NULL
  se.e0<-NULL
  temp<-rep(0,n)
  nb<-rep(0,n)
  e0b<-rep(0,nboot)
  call <- match.call()
  x <- as.matrix(x)
  n <- length(y)
  saveii <- NULL
  fit0 <- theta.fit(x, y, ...)
  yhat0 <- theta.predict(fit0, x, ...)
  app.err <- mean(err.meas(y, yhat0))
  err1 <- matrix(0, nrow = nboot, ncol = n)
  err3 <- err1
  err2 <- rep(0, nboot)

  #compute matrix of bootstrap subsamples
  if(calce00 == F) {
    saveii <- matrix(sample(rep(1:n, nboot),
                            replace = F, size = n * nboot),
                    ncol = nboot, byrow = F)
  }
  # modifies the bootstrap sampling to always use the
  # same number of distinct training points (and test points)
  if(calce00 == T) {
    saveii <- matrix(0, ncol = nboot, nrow = n)
    nn <- trunc(0.6320 * n)
    for(b in 1:nboot) {
      o <- sample(1:n, size = nn, replace = F)
      saveii[, b] <- sample(o, size = n, replace = T)
    }
  }

  # carry out + store bootstrap model fitting and evaluations
  # cat("fitting")
  for(b in 1:nboot) {
    ii <- saveii[, b]
    fit <- theta.fit(x[ii,  ], y[ii], ...)
    yhat1 <- theta.predict(fit, x[ii,  ], ...)
    yhat2 <- theta.predict(fit, x, ...) #what is the "ii"?
    err1[b,  ] <- err.meas(y, yhat2)
    err2[b] <- mean(err.meas(y[ii], yhat1))
    err3[b,  ] <- err.meas(y[ii], yhat1)

    # qhat<-matrix((table(c(y[ii],unique(y)))-1)/length(y), nrow=1)
    # dd<-phat%*%qhat
  }

  #error - error on in-sample points
  #opold <- mean(apply(err1, 1, mean) - err2)

  #u0<-mean(apply(err1, 1, mean))

  #helper to count number of instances of i in vector
  #of indices x
  helper <- function(x, i) {
    sum(x == i)
  }

  #compute OOB error
  #e1<-0
  #dd<-rep(0,n)
  e0i<-rep(0,n) #OOB error for each obs
  for(i in 1:n) {
    o <- apply(saveii, 2, helper, i)
    #nb[i]<-sum(o==0)
    if(sum(o == 0) == 0)
      cat("increase nboot for computation of the .632 estimator",
          fill = T)
    e0i[i]<- sum(err1[o == 0, i])/sum(o == 0)
  }
  e0 <- mean(e0i) #OOB error

  #calcluate SE estimates
  if(calcse){

    count_mat <- matrix(0, n, nboot)
    for(i in 1:n) {
      count_mat[i, ] <- apply(saveii, 2, helper, i)
    }

    #compute terms in (39)
    #cat("39")
    en<-(1-1/n)^(-n)
    qq<-rep(0,nboot)
    nq<-matrix(0,nrow=nboot,ncol=n)

    for(b in 1:nboot){
      for(i in 1:n){
        nq[b,i]<-sum(saveii[,b]==i)
    }}

    for(b in 1:nboot){ #terms in (36)
      oo<-length(unique(saveii[,b]))
      o<-!duplicated(c(unique(saveii[,b]),1:n))
      imiss<-(1:n)[o[-(1:oo)]]
      qq[b]<-sum(err1[b,imiss])/n
    }

    dd<-rep(0,n)
    cov<-function(x,y){cor(x,y)*sqrt(var(x))*sqrt(var(y))}
    for(i in 1:n){
      dd[i]<-(2+1/(n-1))*((e0i[i]-e0)/n)+en*cov(nq[,i],qq)
    }

    #compute (35)
    se.e0<-sqrt(sum(dd^2))

    #compute line (42)
    # for(b in 1:nboot){
    #   for(i in 1:n) {
    #     # o <- apply(saveii, 2, helper, i)
    #     o <- count_mat[i, ]
    #     o[b]<-1
    #     if(sum(o == 0) == 0)
    #       cat("increase nboot for computation of the .632 estimator",
    #           fill = T)
    #     temp[i]<- sum(err1[o == 0, i])/sum(o == 0)
    #     e0b[b] <- e0b[b] + temp[i]/n
    # }}

    #compute line (42)
    #cat("42")
    dib <- matrix(0, nrow = nboot, ncol = n)
    for(b in 1:nboot){
      #analog of e0i, leaving out fold b
      e0i_in<-rep(0,n) #OOB error for each obs, without fold b
      for(i in 1:n) {
        o <- count_mat[i, ]
        o[b]<-1
        #nb[i]<-sum(o==0)
        if(sum(o == 0) == 0)
          cat("increase nboot for computation of the .632 estimator",
              fill = T)
        e0i_in[i]<- sum(err1[o == 0, i])/sum(o == 0)
      }
      e0_in <- mean(e0i_in) #OOB error without fold b

      for(i in 1:n){
        dib[b, i]<-(2+1/(n-1))*((e0i_in[i]-e0_in)/n)+en*cov(nq[-b,i],qq[-b])
      }
    }
    sein.e0 <- sqrt(sum(apply(dib, 2, sd)^2*(nboot-1)^2/nboot))

    #sein.e0<-sqrt(((nboot-1)/nboot)*sum( (e0b-mean(e0b))^2) )

  }

  op632 = 0.632 * (e0 - app.err)
  se.naive=sqrt(var(e0i)/n)
  se.e0.adj = sqrt(se.e0^2 - sein.e0^2) #line (43)
  ratio = (app.err + op632) / e0 #heuristic adjustment to get 632 SE

  return(list("err_hat" = app.err + op632,
              "bias_est" = e0 - app.err - op632,
              "raw_mean" = e0,
              "se_naive" = se.naive,
              "se_est" = se.e0 * ratio,
              "se_est2" = se.e0.adj * ratio,
              "ci_lo" = app.err + op632 - qnorm(1-alpha/2) * se.e0.adj * ratio,
              "ci_hi" = app.err + op632 + qnorm(1-alpha/2) * se.e0.adj * ratio))
}
