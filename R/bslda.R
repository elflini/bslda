bslda <- function(x,y,network=NA,alpha=0.5) {

  covx <- function(x,y, network=NA) {

    nc <- length(unique(y))
    gene.list <- colnames(x)

    if(is.list(network)) {
      prior.list <- unlist(network)
      p.rest <- setdiff(gene.list,prior.list)

      if (length(p.rest)==0){
        path <- network;
      } else {
        path <- c(network,list(p.rest))
      }
    } else {
      path <- list(gene.list)

    }

    L <- length(path)

    if(is.list(network)) {
      scov <- matrix(0,length(gene.list),length(gene.list))
      colnames(scov) <- rownames(scov) <- gene.list
      for (l in 1:L) {
        x.temp <- x[,match(path[[l]],gene.list)]

        ss <- 0
        for (kk in 1:nc) {
          x.temp.y <- x.temp[y==kk,]
          nk <- nrow(x.temp.y)-1
          ss <- ss + cov(x.temp.y)*nk
        }
        scov[match(path[[l]],gene.list),match(path[[l]],gene.list)]<- ss/(length(y)-nc)
      }

    } else {

      for (l in 1:L) {
        x.temp <- x[,match(path[[l]],gene.list)]

        ss <- 0
        for (kk in 1:nc) {
          x.temp.y <- x.temp[y==kk,]
          nk <- nrow(x.temp.y)-1
          ss <- ss + cov(x.temp.y)*nk
        }
        scov <- ss/(length(y)-nc)
      }
    }
    scov+diag(rep(1,ncol(scov)))
  }

  sigma <- covx(x,y,network)

  nc <- length(unique(y))
  L <- length(sigma)

  y1 <- factor(as.numeric(y),levels=1:nc)

  xk <- pk <- muk <- px <- list()
  for (kk in 1: nc) {
    xk[[kk]] <- x[y1==kk,]
    pk[[kk]] <- sum(y1==kk)/length(y1)
    muk[[kk]] <- apply(xk[[kk]],2,mean)
    px[[kk]] <- sweep(x,2,muk[[kk]]/2)
  }

  if (nc==2) {

    mu2 <- (muk[[1]]+muk[[2]])/2
    px2 <- sweep(x,2,mu2)
    mu2r <- muk[[1]]-muk[[2]]

    fit <- enet(sigma,mu2r,lambda=alpha,normalize=F,intercept=F)
    coeff <- predict(fit,s=alpha, type="coefficients",mode="fraction")$coefficients
    pred <- predict(fit,px2,s=alpha, type="fit",mode="fraction")$fit
    dx <- pred+log(pk[[1]]/pk[[2]])

    fitt.v <- factor(ifelse(dx>0,1,2),levels=1:2)
    ac <- sum(diag(table(y1,fitt.v)))/length(y1)

    bac <- NULL
    for (i in 1:nc){
      bac[i] <- sum(y1[y1==i]==fitt.v[y1==i])/length(y1[y1==i])
    }
    bac <- mean(bac)
    #print(ac)

  } else {

    dx <- matrix(0,nrow=nrow(x),ncol=nc)
    coeff <- matrix(0,nrow=ncol(x),ncol=nc)
    rownames(coeff) <- colnames(sigma)

    fit <- pred <- fit.cv <- list()

    for (kk in 1: nc) {

      muk.k <- muk[[kk]]
      px.k <- px[[kk]]

      fit[[kk]] <- enet(sigma,muk.k,lambda=alpha,normalize=F,intercept=F)
      coeff[,kk] <- predict(fit[[kk]],s=alpha, type="coefficients",mode="fraction")$coefficients
      pred[[kk]] <- predict(fit[[kk]],px.k,s=alpha, type="fit",mode="fraction")$fit
      dx[,kk] <- pred[[kk]]+log(pk[[kk]])
    }


    fitt.v <- factor(apply(dx,1,which.max),levels=1:nc)
    ac <- sum(diag(table(y1,fitt.v)))/length(y1)

    bac <- NULL
    for (i in 1:nc){
      bac[i] <- sum(y1[y1==i]==fitt.v[y1==i])/length(y1[y1==i])
    }
    bac <- mean(bac)
  }
  #print(ac)

  re_data <- list(x=x, y=c(y))
  re_estimate <- list(Sigma=sigma, Mu=muk, Pi=pk, alpha=alpha)
  re_result <- list(Fit=fit,pred=pred, coef=coeff, dx=dx, fitted=as.numeric(fitt.v), ac=ac, bac=bac)

  methods::new("bslda",
               data = re_data,
               init = re_estimate,
               result = re_result)

}


