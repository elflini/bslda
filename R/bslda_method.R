
# generic methods for "bslda" class

setMethod(
  f="show",
  signature="bslda",
  definition=function( object ) {

    data <- object@data
    init <- object@init
    result <- object@result

    cat( "Summary: bslda (class: bslda)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Model settings \n")
    cat( "Number of samples to be analyzed: ", nrow(data$x), "\n", sep="" )
    cat( "Number of variables to be analyzed: ", ncol(data$x), "\n", sep="" )
    cat( "Number of classes: ", length(unique(data$y)), "\n", sep="" )
    cat( "Weight of L2 norm: ", init$alpha, "\n", sep="" )
    cat( "Number of non-zero coefficients: ", sum(result$coef!=0), "\n", sep="" )
    cat( "Classification accuracy: ", result$ac, "\n", sep="" )
    cat( "Balanced classification accuracy: ", result$bac, "\n", sep="" )
    cat( "--------------------------------------------------\n" )

  }
)

setGeneric("coeff", def=function(object) standardGeneric("coeff"))

setMethod(
  f="coeff",
  signature="bslda",
  definition=function( object ) {

    # extract objects

    result <- object@result

    coef <- result$coef
    ac <- result$ac
    bac <- result$bac

    return(list(
      coef = coef,
      ac = ac,
      bac=bac

    ))
  }
)



setMethod(
  f="predict",
  signature="bslda",
  definition=function( object,x,y=NA ) {

    # extract objects

    nc <- length(unique(object@data$y))
    L <- length(object@init$Sigma)
    muk <- object@init$Mu
    pi <- object@init$Pi
    discriminant <- object@result$coef

    if (!all(is.na(y))) y1 <- factor(as.numeric(y),levels=1:nc)

    px <- list()
    for (kk in 1: nc) {
      px[[kk]] <- sweep(x,2,muk[[kk]]/2)
    }

    if (nc==2) {
      discriminant <- discriminant[match(colnames(x),names(discriminant))]

      mu2 <- (muk[[1]]+muk[[2]])/2
      px2 <- sweep(x,2,mu2)

      dx <- px2%*%discriminant + log(pi[[1]]/pi[[2]])
      pred <- ifelse(dx>0,1,2)

    } else {
      discriminant <- discriminant[match(colnames(x),rownames(discriminant)),]

      dx <- matrix(,nrow=nrow(x),ncol=nc)
      for (kk in 1:nc) {
        dx[,kk] <- px[[kk]]%*%discriminant[,kk] + log(pi[[kk]])
      }
      pred <- factor(apply(dx,1,which.max),levels=1:nc)
    }

    if (all(is.na(y))) {
      result <- cbind(dx,pred)
    } else {
      ac <- sum(diag(table(y1,pred)))/length(y1)

      bac <- NULL
      for (i in 1:nc){
        bac[i] <- sum(y1[y1==i]==pred[y1==i])/length(y1[y1==i])
      }
      bac <- mean(bac)

      result <- list(pred=cbind(dx,pred), ac=ac, bac=bac)
    }

    return(result)
  }
)

