#' 
#' @export
cv.fGMD <- function(x, y, group,alpha, GGder, lambdader ,lambda = NULL, pred.loss = c("misclass", 
    "loss", "L1", "L2"), nfolds = 5, foldid, delta, lambda.factor,
    nlambda,nfolder,nalpha,nlamder,lamdermin, lamdermax,alphamin, alphamax, loss=loss, ...) {
       if (missing(pred.loss)) 
        pred.loss <- "default" else pred.loss <- match.arg(pred.loss)
    N <- nrow(x)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
    if (missing(delta)) 
        delta <- 1
    if (delta < 0) 
        stop("delta must be non-negtive")
    #########################
    
    
        

    if ( is.null(alpha) | is.null(lambdader) )
    {    
        
      if( !is.null(alpha) | !is.null(lambdader))   {par(mfrow=c(1,1))} else {par(mfrow=c(1,2))}
        
    bn <- as.integer(max(group))
    bs <- as.integer(as.numeric(table(group)))
    ix <- rep(NA, bn)
    iy <- rep(NA, bn)
    j <- 1
    for (g in 1:bn) {
        ix[g] <- j
        iy[g] <- j + bs[g] - 1
        j <- j + bs[g]
    }
    ix <- as.integer(ix)
    iy <- as.integer(iy)
    maximlam=0
    vl=y%*%x;
    for( g in 1:bn )
    {
        tem = vl[ix[g]:iy[g]]
        maximlam = max(MFSGrp::normcpp(tem, diag(1)), maximlam)
    }
    #print(lambda.factor)
    maximlam=maximlam/length(y)
    lam=exp(seq(log(maximlam),log(lambda.factor*maximlam), length.out=10))
    #print(lam)
    
    
    
    #nfolder=3
    foldider <- sample(rep(seq(nfolder), length = N))
    
    #############################
    
    if(is.null(lambdader))
    {
        if (is.null(alpha)) {alp=0;} else {alp=alpha}
        #lam=(1-alp)*lam
        lambdaders=exp(seq(log(lamdermin), log(lamdermax),len=nlamder))
        #lambdaders=append(0, lambdaders)
        #nlamder=nlamder+1
        #nlamder=10
        lambdadersmse=rep(NA, nlamder)
        systimes=rep(NA, nlamder)
    
    
    
    for (j in seq(nlamder)) 
    {          start_time <- Sys.time()
          outlister <- as.list(seq(nfolder))
        for (i in seq(nfolder)) {
            cat("       of lambdas of ", length(lam), " in ",i,"th fold of ", nfolder ," for the ", j, "th lambdader of", nlamder  , "\r" )
            which <- foldider == i
            y_sub <- y[!which]
            outlister[[i]] <- fGMD(x = x[!which, , drop = FALSE],alpha=alp,GGder=GGder,lambdader=lambdaders[j], y = y_sub, group = group, 
                                      lambda = lam, delta = delta,lambda.factor, ...)

            #print(outlister[[i]])
        }
  
        cvstuffer= cv.logit(outlister, lam, x, y, foldider, pred.loss, delta) 
        cvmer <- cvstuffer$cvm
        #print(min(cvstuffer$cvm))
        #cat(j, "\n")
        #cat(lambdaders[j])
        systimes[j]=Sys.time()-start_time
        lambdadersmse[j]=(min(cvmer))
    }

#print(lambdadersmse)
#print(lambdaders)
#print(systimes)
       # lamders[1]=exp(-50)
plot(y=lambdadersmse, x=log(lambdaders )) 
chose= which(lambdadersmse==min(lambdadersmse) )[1]
#print(chose[1])
lambdader=min(lambdaders[chose])
#lambdader=5e-4
#cat(systimes[chose[1]])
ET=nfolds*nlambda*systimes[chose]/(nfolder*nlamder)
#cat(systimes[chose], "\n")
cat("\r        Chosen lambdader is", lambdader, "and Maximum Estimated Time:",  ET  ," seconds           \r \n" )
#stop("here")
    
    }
   ############################### 
    
   
    if(is.null(alpha))
    {
        
        #nalpha=9
        alphasmse=rep(NA, nalpha)
        nalphaseq=nalpha+2
        alphas=seq(alphamin, alphamax, len=nalphaseq )[-c(1,nalphaseq)]
        #alphas=seq(0.1,0.9, len=nalpha )
        systimes=rep(NA, nalpha)
        
        
        for (j in seq(nalpha)) 
        {          start_time <- Sys.time()
        outlister <- as.list(seq(nfolder))
        for (i in seq(nfolder)) {
            cat("       of lambdas of ", length(lam), " in ",i,"th fold of ", nfolder ," for the ", j, "th alpha of", nalpha  , "\r" )
            which <- foldider == i
            y_sub <- y[!which]
            outlister[[i]] <- fGMD(x = x[!which, , drop = FALSE],alpha=alphas[j],GGder=GGder,lambdader=lambdader, y = y_sub, group = group, 
                                      lambda = lam, delta = delta,lambda.factor, ...)
            
            #print(outlister[[i]])
        }
        
        cvstuffer= cv.logit(outlister, lam, x, y, foldider, pred.loss, delta) 
        cvmer <- cvstuffer$cvm
        #print(min(cvstuffer$cvm))
        #cat(j, "\n")
        #cat(lambdaders[j])
        systimes[j]=Sys.time()-start_time
        alphasmse[j]=(min(cvmer))
        }
        
        #print(lambdadersmse)
        #print(lambdaders)
        #print(systimes)
        
        plot(y=alphasmse, x=alphas )
        chose= which(alphasmse==min(alphasmse) )[1]
        #print(chose[1])
        alpha=min(alphas[chose])
        #lambdader=5e-4
        #cat(systimes[chose[1]])
        ET=nfolds*nlambda*systimes[chose]/(nfolder*nalpha)
        #cat(systimes[chose], "\n")
        cat("\r        Chosen alpha is", alpha, "and Maximum Estimated Time:",  ET  ," seconds                    \r \n" )
        #stop("here")
        
    }
    ############################### 
    
    
    
     
    
    
    }  
    
    
    
    
    
    
    
    
    
    ######################
    
    fGMD.object <- fGMD(x, y, group,alpha=alpha, GGder=GGder,lambdader=lambdader,lambda = lambda, delta = delta, lambda.factor=lambda.factor,nlambda=nlambda, ...)
    lambda <- fGMD.object$lambda
    #print(lambda)
    # predict -> coef
    if (missing(foldid)) 
        foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist <- as.list(seq(nfolds))
    ###Now fit the nfold models and store them
    for (i in seq(nfolds)) {
        cat("       of lambdas out of ", length(lambda), " grid points in the ", i,"th fold of ",nfolds," folds. ",  "\r")
        which <- foldid == i
        y_sub <- y[!which]
        outlist[[i]] <- fGMD(x = x[!which, , drop = FALSE],alpha=alpha,GGder=GGder,lambdader=lambdader, y = y_sub, group = group, 
            lambda = lambda, delta = delta,lambda.factor=lambda.factor,nlambda=nlambda, ...)
    }
    ###What to do depends on the pred.loss and the model fit
    fun <- paste("cv", class(fGMD.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, pred.loss, delta))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, 
        cvlo = cvm - cvsd, name = cvname, fGMD.fit = fGMD.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.fGMD"
    obj
}


#' @export
cv.hsvm <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    if (pred.loss == "default") 
        pred.loss <- "misclass"
    if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
        warning("Only 'misclass' and 'loss' available for HHSVM classification; 'misclass' used")
        pred.loss <- "misclass"
    }
    ###Turn y into c(0,1)
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- switch(pred.loss, loss = 2 * hubercls(y * predmat, delta), misclass = (y != 
        ifelse(predmat > 0, 1, -1)))
	N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}



#' @export
cv.logit <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    if (pred.loss == "default") 
        pred.loss <- "misclass"
    if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
        warning("Only 'misclass' and 'loss' available for logistic regression; 'misclass' used")
        pred.loss <- "misclass"
    }
    prob_min <- 1e-05
    fmax <- log(1/prob_min - 1)
    fmin <- -fmax
    ###Turn y into c(0,1)
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    predmat <- pmin(pmax(predmat, fmin), fmax)
    cvraw <- switch(pred.loss, loss = 2 * log(1 + exp(-y * predmat)), misclass = (y != 
        ifelse(predmat > 0, 1, -1)))
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}

#' @export
cv.sqsvm <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    if (pred.loss == "default") 
        pred.loss <- "misclass"
    if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
        warning("Only 'misclass' and 'loss' available for squared SVM classification; 'misclass' used")
        pred.loss <- "misclass"
    }
    ###Turn y into c(0,1)
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- switch(pred.loss, loss = 2 * ifelse((1 - y * predmat) <= 0, 0, 
        (1 - y * predmat))^2, misclass = (y != ifelse(predmat > 0, 1, -1)))
	N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}

#' @export
cv.ls <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(L2 = "Least-Squared loss", L1 = "Absolute loss")
    if (pred.loss == "default") 
        pred.loss <- "L2"
    if (!match(pred.loss, c("L2", "L1"), FALSE)) {
        warning("Only 'L2' and 'L1'  available for least squares models; 'L2' used")
        pred.loss <- "L2"
    }
    predmat <- matrix(NA, length(y), length(lambda))
    nfolds <- max(foldid)
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE])
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- switch(pred.loss, L2 = (y - predmat)^2, L1 = abs(y - predmat))
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 
