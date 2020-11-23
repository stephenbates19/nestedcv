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

  if(calce00 == F) {
    saveii <- matrix(sample(rep(1:n, nboot), replace = F, size = n *
                              nboot), ncol = nboot, byrow = F)
  }
  if(calce00 == T) {
    saveii <- matrix(0, ncol = nboot, nrow = n)
    nn <- trunc(0.6320 * n)
    for(b in 1:nboot) {
      o <- sample(1:n, size = nn, replace = F)
      saveii[, b] <- sample(o, size = n, replace = T)
    }
  }
  #
  for(b in 1:nboot) {
    ii <- saveii[, b]
    fit <- theta.fit(x[ii,  ], y[ii], ...)
    yhat1 <- theta.predict(fit, x[ii,  ], ...)
    yhat2 <- theta.predict(fit, x, ...) #what is the "ii"?
    err1[b,  ] <- err.meas(y, yhat2)
    err2[b] <- mean(err.meas(y[ii], yhat1))
    err3[b,  ] <- err.meas(y[ii], yhat1)
    qhat<-matrix((table(c(y[ii],unique(y)))-1)/length(y), nrow=1)
    dd<-phat%*%qhat

  }


  opold <- mean(apply(err1, 1, mean) - err2)

  u0<-mean(apply(err1, 1, mean))
  junk <- function(x, i)
  {
    sum(x == i)
  }
  e0 <- 0
  e1<-0


  dd<-rep(0,n)
  e0i<-rep(0,n)
  for(i in 1:n) {
    o <- apply(saveii, 2, junk, i)
    nb[i]<-sum(o==0)
    if(sum(o == 0) == 0)
      cat("increase nboot for computation of the .632 estimator",
          fill = T)
    e0i[i]<- sum(err1[o == 0, i])/sum(o == 0)
    e0 <- e0 + e0i[i]/n
  }

  if(calcse){

    count_mat <- matrix(0, n, nboot)
    for(i in 1:n) {
      count_mat[i, ] <- apply(saveii, 2, junk, i)
    }

    for(b in 1:nboot){
      for(i in 1:n) {
        # o <- apply(saveii, 2, junk, i)
        # o[b]<-1
        o <- count_mat[, i]
        if(sum(o == 0) == 0)
          cat("increase nboot for computation of the .632 estimator",
              fill = T)
        temp[i]<- sum(err1[o == 0, i])/sum(o == 0)
        e0b[b] <- e0b[b] + temp[i]/n
      }}

    sein.e0<-sqrt(((nboot-1)/nboot)*sum( (e0b-mean(e0b))^2) )



    en<-(1-1/n)^(-n)
    qq<-rep(0,nboot)
    nq<-matrix(0,nrow=nboot,ncol=n)

    for(b in 1:nboot){
      for(i in 1:n){
        nq[b,i]<-sum(saveii[,b]==i)
      }}

    for(b in 1:nboot){
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
    se.e0<-sqrt(sum(dd^2))
  }

  op632 = 0.632 * (e0 - app.err)
  se.naive=sqrt(var(e0i)/n)

  return(list("err_hat" = app.err + op632,
              "bias_est" = e0 - app.err - op632,
              "raw_mean" = e0,
              "se_naive" = se.naive,
              "se_est" = se.e0,
              "ci_lo" = app.err + op632 - qnorm(1-alpha/2) * se.e0,
              "ci_hi" = app.err + op632 + qnorm(1-alpha/2) * se.e0))
}
